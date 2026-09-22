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

N = 100
r_min, r_max = 0.5, 6.0
r = np.linspace(r_min, r_max, N)
dr = r[1] - r[0]
G_target = np.array([Gamma_fib(rr) for rr in r])

# Precompute a_over_nu at each r
ani_arr = np.array([a_over_nu(rr) for rr in r])

def RHS(G):
    out = np.zeros(N)
    # Vectorized
    Gp = np.gradient(G, dr)
    Gpp = np.gradient(Gp, dr)
    for i in range(1, N-1):
        out[i] = (Gpp[i]/r[i] - Gp[i]/(r[i]*r[i])) + ani_arr[i]*Gp[i]
    return out

def evolve(G, dt=1e-4, steps=5000):
    for s in range(steps):
        k1 = RHS(G); k2 = RHS(G + 0.5*dt*k1)
        G = G + dt*k2
        G[0] = G_target[0]; G[-1] = G_target[-1]
    return G

print("="*60)
print("PROBE #2 (fast): dynamics stability")
print("="*60)
print()
print(f"{'mode':>5} {'amp':>7} {'err_init':>14} {'err_final':>14} {'ratio':>10}")
print("-"*60)

for mode in [1, 2]:
    for amp in [0.01, 0.1]:
        pert = amp*np.sin(mode*math.pi*(r-r_min)/(r_max-r_min))
        G_init = G_target + pert
        G_final = evolve(G_init.copy(), steps=5000)
        e0 = np.max(np.abs(G_init - G_target))
        e1 = np.max(np.abs(G_final - G_target))
        ratio = e1/e0 if e0 > 0 else 0
        print(f"{mode:>5} {amp:>7.3f} {e0:>14.4e} {e1:>14.4e} {ratio:>10.4f}")

print()
print("ratio < 1: STABLE (perturbation decays)")
print("ratio > 1: UNSTABLE (grows)")
print("ratio ~ 1: MARGINAL")
