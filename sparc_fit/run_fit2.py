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

# === Fibonacci halo as DENSITY ===
# rho(r) = rho_0 * Gamma_fib(r)
# Gamma_fib(r) = sum_k F_k * exp(-k*lambda*r/rc), lambda = log(phi)
# M(r) = 4*pi*rho_0 * integral_0^r Gamma(r') r'^2 dr'
# v^2(r) = G*M(r)/r

def mass_fib(r, rc):
    """Vectorized M(r) for Fibonacci density profile."""
    r = np.asarray(r)
    lam = LAM / rc  # lambda with rc
    M = np.zeros_like(r)
    for k in range(K):
        if k == 0:
            # integral of r^2 dr = r^3/3
            M += FIB[k] * r**3 / 3.0
        else:
            kl = k * lam
            # integral = 2/kl^3 - exp(-kl r) * (r^2/kl + 2r/kl^2 + 2/kl^3)
            term = 2.0/kl**3 - np.exp(-kl*r)*(r**2/kl + 2*r/kl**2 + 2/kl**3)
            M += FIB[k] * term
    return 4.0 * math.pi * M

def vdm2_fib(data, M0, rc):
    """v^2_DM = G * M(r) / r, M(r) = M0 * (integral part)"""
    r = data[:,0]
    M_unit = mass_fib(r, rc)  # in units of rho_0
    # M0 is total mass scale (rho_0 * 4pi absorbed)
    M = M0 * M_unit / (4.0 * math.pi)  # remove 4pi, M0 = rho_0
    return G_NEWTON * M / np.maximum(r, 1e-9)

# Alternative: treat M_tot directly
def vdm2_fib_v2(data, M_tot, rc):
    """Same but M_tot = total mass at r -> inf"""
    r = data[:,0]
    # Asymptotic mass: integral_0^inf Gamma(r) r^2 dr
    # for k=0: divergent! No - wait, sum F_k exp(-k lam r), F_0=1 (k=0 constant term)
    # So integral_0^inf r^2 dr diverges. Need cutoff. Hmm.
    # Actually k=0 term is constant: rho -> rho_0 as r->0? Wait, rho(0) = sum F_k = 12.
    # Then rho(r) decays as r->inf. But k=0 term gives exp(0) = 1 constant forever!
    # So density -> rho_0*F_0 = rho_0 as r->inf. Total mass = divergent.
    # This means we need K with F_0 removed OR treat F_0 as constant background.
    # Simplest fix: shift so density goes to 0. Use only k>=1 terms in halo.
    r_safe = np.maximum(r, 1e-9)
    lam = LAM / rc
    M_unit = np.zeros_like(r)
    for k in range(1, K):  # skip k=0 (constant term -> divergent mass)
        kl = k * lam
        term = 2.0/kl**3 - np.exp(-kl*r_safe)*(r_safe**2/kl + 2*r_safe/kl**2 + 2/kl**3)
        M_unit += FIB[k] * term
    return G_NEWTON * M_tot * M_unit / r_safe

def vdm2_vortex(data, A, rc):
    x = np.maximum(data[:,0]/rc, 1e-9)
    return A*A*(x - np.arctan(x))/x

def nfw_vsq(data, M200, c):
    r = data[:,0]
    rho_c = 2.775e2
    r200 = (3*M200/(4*math.pi*200*rho_c))**(1/3)
    rs = r200 / max(c, 0.1)
    x = np.maximum(r/rs, 1e-6)
    Mr = M200*(np.log(1+x) - x/(1+x))/(np.log(1+c) - c/(1+c))
    return G_NEWTON * Mr / r

def chi2_single(data, A, rc, ups, model):
    if A <= 0 or rc <= 0 or ups <= 0:
        return 1e30
    vbar = vbar2(data, ups)
    try:
        if model == 'vortex':
            vdm = vdm2_vortex(data, A, rc)
        elif model == 'fib':
            vdm = vdm2_fib_v2(data, A*1e11, rc)
        elif model == 'nfw':
            vdm = nfw_vsq(data, A*1e10, rc)
        vmod = np.sqrt(np.maximum(vbar + vdm, 0))
        res = (data[:,1] - vmod) / np.maximum(data[:,2], 0.01)
        return float(np.sum(res**2))
    except Exception:
        return 1e30

def fit_one(data, model):
    A_grid = np.logspace(0.7, 2.6, 20)
    rc_grid = np.logspace(-0.7, 1.5, 20)
    ups_grid = [0.3, 0.5, 0.8, 1.0]
    best = (1e30, None)
    for A in A_grid:
        for rc in rc_grid:
            for ups in ups_grid:
                c2 = chi2_single(data, A, rc, ups, model)
                if c2 < best[0]:
                    best = (c2, (A, rc, ups))
    if best[1] is not None:
        A0, rc0, u0 = best[1]
        for it in range(3):
            s = 0.3 ** (it+1)
            for dA in [-s*A0, 0, s*A0]:
                for drc in [-s*rc0, 0, s*rc0]:
                    for du in [-0.2*s, 0, 0.2*s]:
                        A = max(A0 + dA, 1e-3)
                        rc = max(rc0 + drc, 1e-3)
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
print()
print(f"{'galaxy':<12} {'N':>4} {'c2VORTEX':>10} {'c2FIB':>10} {'c2NFW':>10}")
print("-" * 60)

sum_v = 0; sum_f = 0; sum_n = 0; n_dof = 0
for name in test_names:
    data = gals[name]
    n = len(data)
    if n < 5:
        continue
    c2_v, _ = fit_one(data, 'vortex')
    c2_f, _ = fit_one(data, 'fib')
    c2_n, _ = fit_one(data, 'nfw')
    dof = max(n - 3, 1)
    print(f"{name:<12} {n:>4} {c2_v/dof:>10.4f} {c2_f/dof:>10.4f} {c2_n/dof:>10.4f}")
    sum_v += c2_v; sum_f += c2_f; sum_n += c2_n; n_dof += dof

print("-" * 60)
print(f"{'TOTAL mean':<12} {'':>4} {sum_v/n_dof:>10.4f} {sum_f/n_dof:>10.4f} {sum_n/n_dof:>10.4f}")
