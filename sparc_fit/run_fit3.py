import numpy as np
import math, os

PHI = (1+math.sqrt(5))/2
LAM = math.log(PHI)
FIB = [1,1,2,3,5,8,13,21,34,55,89,144]
K = 5
G_NEWTON = 4.30091e-6

def parse_massmodels(path):
    galaxies = {}
    with open(path, 'r') as f:
        for line in f:
            if not line or line[0] == '#' or len(line) < 40:
                continue
            try:
                gid = line[0:11].strip()
                if not gid or gid.startswith('Title') or gid.startswith('Author'):
                    continue
                D = float(line[12:18].strip())
                R = float(line[19:25].strip())
                Vobs = float(line[26:32].strip())
                eV = float(line[33:38].strip())
                Vg = float(line[39:45].strip())
                Vd = float(line[46:52].strip())
                Vb = float(line[53:59].strip())
            except (ValueError, IndexError):
                continue
            galaxies.setdefault(gid, []).append((R, Vobs, eV, Vg, Vd, Vb))
    return {k: np.array(v) for k, v in galaxies.items()}

def vbar2(data, ups):
    return data[:,3]**2 + ups*data[:,4]**2 + 1.4*ups*data[:,5]**2

# === GravVortex parametrization (from GravVortexIntegration-v3) ===
# Gamma(r) = Gamma_0 * sum F_k exp(-k*LAM*r/rc)
# omega(r) = |dGamma/dr| / (2*pi*r)
# rho(r) = kappa * omega^2
# M(r) = integral 4*pi*rho*r^2 dr (numeric)
# v_DM^2 = G*M(r)/r

def Gamma_fib(r, rc):
    return sum(FIB[k]*np.exp(-k*LAM*r/rc) for k in range(K))

def dGamma_fib(r, rc, h=1e-6):
    return (Gamma_fib(r+h, rc) - Gamma_fib(r-h, rc))/(2*h)

def rho_omega(r_arr, rc, kappa):
    """rho(r) = kappa * omega^2, omega = |dGamma/dr|/(2*pi*r)"""
    rho = np.zeros_like(r_arr)
    for i, r in enumerate(r_arr):
        if r < 1e-9:
            continue
        dg = dGamma_fib(r, rc)
        om = abs(dg) / (2*math.pi*r)
        rho[i] = kappa * om*om
    return rho

def M_rho(r_max, rc, kappa, n_int=200):
    """M(r) = integral_0^r 4*pi*rho*r^2 dr via trapezoid"""
    r_int = np.linspace(1e-4, r_max, n_int)
    rho = rho_omega(r_int, rc, kappa)
    integrand = 4*math.pi*rho*r_int**2
    M = np.trapz(integrand, r_int)
    return M

def vdm2_gravvortex(data, A_scale, rc):
    """A_scale absorbs Gamma_0^2 * kappa. v^2 = G*M/r"""
    r = data[:,0]
    vdm = np.zeros_like(r)
    for i, ri in enumerate(r):
        if ri < 1e-9:
            continue
        M = M_rho(ri, rc, A_scale, n_int=100)
        vdm[i] = G_NEWTON * M / ri
    return vdm

# For speed: precompute M(r) on a fine grid, then interpolate
def M_table(rc, kappa, r_max=60.0, n_pts=400):
    r_grid = np.linspace(1e-4, r_max, n_pts)
    rho = rho_omega(r_grid, rc, kappa)
    integrand = 4*math.pi*rho*r_grid**2
    M = np.zeros(n_pts)
    for i in range(1, n_pts):
        M[i] = np.trapz(integrand[:i+1], r_grid[:i+1])
    return r_grid, M

def vdm2_gravvortex_fast(data, A_scale, rc, r_grid, M_grid):
    r = data[:,0]
    M_interp = np.interp(r, r_grid, M_grid, left=0.0, right=M_grid[-1])
    return G_NEWTON * (A_scale * M_interp) / np.maximum(r, 1e-9)

def vdm2_vortex(data, A, rc):
    x = np.maximum(data[:,0]/rc, 1e-9)
    return A*A*(x - np.arctan(x))/x

def chi2_single(data, A, rc, ups, model, cache=None):
    if A <= 0 or rc <= 0 or ups <= 0:
        return 1e30
    vbar = vbar2(data, ups)
    try:
        if model == 'vortex':
            vdm = vdm2_vortex(data, A, rc)
        elif model == 'gravvortex':
            r_grid, M_grid = M_table(rc, 1.0)
            # A in grid is Gamma_0^2 * kappa
            vdm = vdm2_gravvortex_fast(data, A, rc, r_grid, M_grid)
        vmod = np.sqrt(np.maximum(vbar + vdm, 0))
        res = (data[:,1] - vmod) / np.maximum(data[:,2], 0.01)
        return float(np.sum(res**2))
    except Exception as e:
        return 1e30

def fit_one(data, model):
    A_grid = np.logspace(0.5, 5.0, 25)  # wider for gravvortex
    rc_grid = np.logspace(-0.7, 1.5, 20)
    ups_grid = [0.3, 0.5, 0.8, 1.0]
    best = (1e30, None)
    for A in A_grid:
        for rc in rc_grid:
            for ups in ups_grid:
                c2 = chi2_single(data, A, rc, ups, model)
                if c2 < best[0]:
                    best = (c2, (A, rc, ups))
    # Local refine
    if best[1] is not None:
        A0, rc0, u0 = best[1]
        for it in range(2):
            s = 0.3 ** (it+1)
            for dA in [-s*A0, 0, s*A0]:
                for drc in [-s*rc0, 0, s*rc0]:
                    for du in [-0.2*s, 0, 0.2*s]:
                        A = max(A0 + dA, 1e-6)
                        rc = max(rc0 + drc, 1e-6)
                        ups = max(u0 + du, 1e-3)
                        c2 = chi2_single(data, A, rc, ups, model)
                        if c2 < best[0]:
                            best = (c2, (A, rc, ups))
                            A0, rc0, u0 = A, rc, ups
    return best

# === Main ===
path = os.path.expanduser('~/sparc-vortex-nfw/data/MassModels_Lelli2016c.mrt')
gals = parse_massmodels(path)
print(f"Loaded {len(gals)} galaxies")
print()

test_names = ['D564-8', 'DDO064', 'F561-1', 'F563-V1', 'F565-V2', 'D631-7', 'DDO170', 'UGCA444']
test_names = [n for n in test_names if n in gals]
print(f"Testing on {len(test_names)} galaxies")
print(f"Model: rho(r) = kappa * omega^2, omega = |dGamma/dr|/(2*pi*r)")
print(f"  (from GravVortexIntegration-v3)")
print()
print(f"{'galaxy':<12} {'N':>4} {'c2VORTEX':>10} {'c2GravVtx':>12} {'c2NFW':>10}")
print("-" * 60)

sum_v = 0; sum_g = 0; n_dof = 0
for name in test_names:
    data = gals[name]
    n = len(data)
    if n < 5:
        continue
    c2_v, _ = fit_one(data, 'vortex')
    c2_g, _ = fit_one(data, 'gravvortex')
    dof = max(n - 3, 1)
    print(f"{name:<12} {n:>4} {c2_v/dof:>10.4f} {c2_g/dof:>12.4f}")
    sum_v += c2_v; sum_g += c2_g; n_dof += dof

print("-" * 60)
print(f"{'TOTAL mean':<12} {'':>4} {sum_v/n_dof:>10.4f} {sum_g/n_dof:>12.4f}")
print()
print("GravVortex: rho = kappa*(|dGamma/dr|/(2*pi*r))^2")
print("VORTEX: standard SPARC reference")
