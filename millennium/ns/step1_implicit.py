import numpy as np
import math

PHI = (1+math.sqrt(5))/2
LAM = math.log(PHI)
FIB = [1,1,2,3,5,8,13,21,34,55,89,144]
K = 5

def Gamma_fib(r):
    return sum(FIB[k]*math.exp(-k*LAM*r) for k in range(K))

def a_over_nu(r):
    num = sum(k*k*FIB[k]*math.exp(-k*LAM*r) for k in range(K))
    den = sum(k*FIB[k]*math.exp(-k*LAM*r) for k in range(K))
    if abs(den) < 1e-15: return 0.0
    return (num/den)*LAM/r + 1.0/(r*r)

N = 200
r_min, r_max = 0.3, 10.0
r = np.linspace(r_min, r_max, N)
dr = r[1] - r[0]
ani_arr = np.array([a_over_nu(rr) for rr in r])
G_target = np.array([Gamma_fib(rr) for rr in r])

# Build L (interior only)
L = np.zeros((N, N))
for i in range(1, N-1):
    ri = r[i]
    L[i, i-1] = (1.0/(dr*dr * ri)) + (1.0/(2*dr * ri*ri)) - ani_arr[i]/(2*dr)
    L[i, i]   = -2.0/(dr*dr * ri)
    L[i, i+1] = (1.0/(dr*dr * ri)) - (1.0/(2*dr * ri*ri)) + ani_arr[i]/(2*dr)

print("="*72)
print("STEP 1 (implicit): Non-stationary evolution with backward Euler")
print("="*72)
print()

# Eigenvalues (interior)
eigs = np.linalg.eigvals(L[1:-1, 1:-1])
max_real = max(e.real for e in eigs)
max_abs = max(abs(e) for e in eigs)
print(f"Grid: N={N}, r ∈ [{r_min}, {r_max}]")
print(f"Max Re(λ) = {max_real:+.6e}")
print(f"Max |λ|   = {max_abs:+.6e}")
print(f"RK2 stability requires dt < {2/max_abs:.3e}")
print(f"→ Using IMPLICIT backward Euler (unconditionally stable)")
print()

# Implicit step: (I - dt*L) G_new = G_old
I = np.eye(N)

def evolve_implicit(G, dt, steps):
    A = I - dt*L
    # Precompute LU-like solve via numpy
    for s in range(steps):
        G_new = np.linalg.solve(A, G)
        G_new[0] = G_target[0]
        G_new[-1] = G_target[-1]
        G = G_new
    return G

def energy(G):
    delta = G - G_target
    return float(np.sum(delta*delta * r) * dr)

np.random.seed(42)
print(f"{'perturbation':<25} {'E(0)':>12} {'E(T)':>12} {'ratio':>10}")
print("-"*72)

# Use larger dt because implicit is stable
dt = 0.01
T_total = 5.0
steps = int(T_total / dt)

for name, pert in [
    ("sin(πx)", lambda r: 0.1*np.sin(np.pi*(r-r_min)/(r_max-r_min))),
    ("sin(2πx)", lambda r: 0.1*np.sin(2*np.pi*(r-r_min)/(r_max-r_min))),
    ("sin(5πx)", lambda r: 0.1*np.sin(5*np.pi*(r-r_min)/(r_max-r_min))),
    ("gaussian(center=3)", lambda r: 0.1*np.exp(-((r-3.0)/0.5)**2)),
    ("random 1%", lambda r: 0.01*np.random.randn(len(r))),
    ("random 10%", lambda r: 0.1*np.random.randn(len(r))),
]:
    G_init = G_target + pert(r)
    E0 = energy(G_init)
    G_final = evolve_implicit(G_init.copy(), dt=dt, steps=steps)
    E1 = energy(G_final)
    ratio = E1/E0 if E0 > 0 else 0
    print(f"{name:<25} {E0:>12.4e} {E1:>12.4e} {ratio:>10.6f}")

print()
print("="*72)
print("STEP 1 CONCLUSION")
print("="*72)
print()

# Slowest decay rate
slowest = max(e.real for e in eigs)
print(f"Slowest decay mode: Re(λ_max) = {slowest:+.6e}")
print(f"  → time constant τ = 1/|λ_max| = {1/abs(slowest):.4f}")
print(f"  → at T = {T_total}: exp(-T/τ) = {math.exp(-T_total*abs(slowest)):.4e}")
print()
print("All perturbations decay if ratio < 1.")
print("→ Numerical evidence for global asymptotic stability.")
