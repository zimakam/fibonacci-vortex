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

def gamma_fib_arr(r, rc):
    return sum(FIB[k]*np.exp(-k*LAM*r/rc) for k in range(K))

def vdm2_vortex(data, A, rc):
    x = np.maximum(data[:,0]/rc, 1e-9)
    return A*A*(x - np.arctan(x))/x

def vdm2_fib(data, M_tot, rc):
    r = data[:,0]
    G0 = sum(FIB[:K])
    G_inf = FIB[0]
    Gf = gamma_fib_arr(r, rc)
    f = (G0 - Gf) / (G0 - G_inf)
    M = M_tot * np.clip(f, 0.0, 1.0)
    return G_NEWTON * M / r

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
    if model == 'vortex':
        vdm = vdm2_vortex(data, A, rc)
    elif model == 'fib':
        vdm = vdm2_fib(data, A*1e10, rc)
    elif model == 'nfw':
        vdm = nfw_vsq(data, A*1e10, rc)
    vmod = np.sqrt(np.maximum(vbar + vdm, 0))
    res = (data[:,1] - vmod) / np.maximum(data[:,2], 0.01)
    return float(np.sum(res**2))

def fit_one(data, model, A_grid=None, rc_grid=None, ups_grid=None):
    if A_grid is None:
        A_grid = np.logspace(0.7, 2.6, 20)
    if rc_grid is None:
        rc_grid = np.logspace(-0.7, 1.5, 20)
    if ups_grid is None:
        ups_grid = [0.3, 0.5, 0.8, 1.0]
    best = (1e30, None)
    for A in A_grid:
        for rc in rc_grid:
            for ups in ups_grid:
                c2 = chi2_single(data, A, rc, ups, model)
                if c2 < best[0]:
                    best = (c2, (A, rc, ups))
    # Local refine: 3 iterations halving step
    if best[1] is not None:
        A0, rc0, u0 = best[1]
        for it in range(3):
            scale = 0.3 ** (it+1)
            for dA in [-scale*A0, 0, scale*A0]:
                for drc in [-scale*rc0, 0, scale*rc0]:
                    for du in [-0.2*scale, 0, 0.2*scale]:
                        A = A0*(1+dA/A0) if A0 > 0 else A0
                        rc = rc0 + drc
                        ups = u0 + du
                        if A <= 0 or rc <= 0 or ups <= 0:
                            continue
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
print(f"{'galaxy':<12} {'N':>4} {'c2VORTEX':>10} {'c2FIB':>10} {'c2NFW':>10}   (chi2/dof)")
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
print()
print("Lower chi2/dof = better fit")
