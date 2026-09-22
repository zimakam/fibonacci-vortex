import numpy as np
import math

PHI = (1+math.sqrt(5))/2
LAM = math.log(PHI)
FIB = [1,1,2,3,5,8,13,21,34,55,89,144]
K = 5

def Gamma_fib(r):
    return sum(FIB[k]*np.exp(-k*LAM*r) for k in range(K))

def Gamma_p(r, h=1e-5):
    return (Gamma_fib(r+h) - Gamma_fib(r-h))/(2*h)

def Gamma_pp(r, h=1e-5):
    return (Gamma_fib(r+h) - 2*Gamma_fib(r) + Gamma_fib(r-h))/(h*h)

def a_over_nu(r):
    """Из нашей identity: a/nu = <k>_w * lambda / r + 1/r^2"""
    num = sum(k*k*FIB[k]*math.exp(-k*LAM*r) for k in range(K))
    den = sum(k*FIB[k]*math.exp(-k*LAM*r) for k in range(K))
    if abs(den) < 1e-15:
        return 0.0
    kbar = num/den
    return kbar*LAM/r + 1.0/(r*r)

# === Numerical evolution ===
# Grid: r in [0.3, 8.0], N points
N = 400
r_min, r_max = 0.3, 8.0
r = np.linspace(r_min, r_max, N)
dr = r[1] - r[0]

# Initial condition: Gamma_fib + small perturbation
def make_ic(perturb_amp=0.01, perturb_mode=1):
    G0 = np.array([Gamma_fib(rr) for rr in r])
    # Perturbation: sin(pi * n * (r-r_min)/(r_max-r_min))
    pert = perturb_amp * np.sin(perturb_mode * math.pi * (r - r_min)/(r_max - r_min))
    return G0 + pert

# Laplacian-like operator from Burgers
def dG_dr(G, i):
    if i == 0:
        return (G[1] - G[0]) / dr
    elif i == N-1:
        return (G[N-1] - G[N-2]) / dr
    else:
        return (G[i+1] - G[i-1]) / (2*dr)

def d2G_dr2(G, i):
    if i == 0:
        return (G[2] - 2*G[1] + G[0]) / (dr*dr)
    elif i == N-1:
        return (G[N-1] - 2*G[N-2] + G[N-3]) / (dr*dr)
    else:
        return (G[i+1] - 2*G[i] + G[i-1]) / (dr*dr)

def RHS(G):
    """RHS = nu*(G''/r - G'/r^2) + a(r)*G'"""
    nu = 1.0  # absorb into a/nu ratio
    out = np.zeros(N)
    for i in range(1, N-1):
        Gp = dG_dr(G, i)
        Gpp = d2G_dr2(G, i)
        ani = a_over_nu(r[i])
        out[i] = nu*(Gpp/r[i] - Gp/(r[i]*r[i])) + nu*ani*Gp
    return out

# Run RK2 evolution
def evolve(G, dt=1e-4, steps=20000):
    history = []
    for s in range(steps):
        k1 = RHS(G)
        k2 = RHS(G + 0.5*dt*k1)
        G_new = G + dt*k2
        # Boundary conditions: fixed at Gamma_fib
        G_new[0] = Gamma_fib(r[0])
        G_new[-1] = Gamma_fib(r[-1])
        G = G_new
        if (s+1) % 2000 == 0:
            err = np.max(np.abs(G - np.array([Gamma_fib(rr) for rr in r])))
            history.append((s*dt, err))
    return G, history

# === Test stability for different perturbations ===
print("="*72)
print("PROBE #2: Dynamics of Gamma_fib under Burgers evolution")
print("="*72)
print()
print("Evolution: d_t Gamma = nu*(G''/r - G'/r^2) + nu*(a/nu)(r)*G'")
print("Stationary solution: Gamma_fib")
print()
print("Question: if we perturb Gamma_fib, does it return?")
print()

print(f"{'mode':>6} {'amp':>8} {'t=0.5':>12} {'t=1.0':>12} {'t=2.0':>12} {'final':>12}")
print("-"*72)

for mode in [1, 2, 3]:
    for amp in [0.01, 0.1]:
        G_init = make_ic(perturb_amp=amp, perturb_mode=mode)
        G_final, hist = evolve(G_init, dt=1e-4, steps=20000)
        err0 = hist[0][1] if hist else 0
        err_mid = hist[2][1] if len(hist) > 2 else 0
        err_1 = hist[4][1] if len(hist) > 4 else 0
        err_final = hist[-1][1] if hist else 0
        print(f"{mode:>6} {amp:>8.3f} {err0:>12.4e} {err_mid:>12.4e} {err_1:>12.4e} {err_final:>12.4e}")

print()
print("Интерпретация:")
print("  Если err -> 0: решение УСТОЙЧИВО (притягивает)")
print("  Если err растёт: НЕУСТОЙЧИВО")
print("  Если err постоянна: МАРГИНАЛЬНО (нейтрально)")

# === Test: evolution starting far from Gamma_fib ===
print()
print("="*72)
print("Test: start from random perturbation of size 1.0")
print("="*72)
np.random.seed(42)
G_init = np.array([Gamma_fib(rr) for rr in r]) + 1.0*np.random.randn(N)*0.1
G_final, hist = evolve(G_init, dt=1e-4, steps=20000)
print(f"Initial error: {np.max(np.abs(G_init - np.array([Gamma_fib(rr) for rr in r]))):.4e}")
print(f"Final   error: {np.max(np.abs(G_final - np.array([Gamma_fib(rr) for rr in r]))):.4e}")
print(f"Ratio: {np.max(np.abs(G_final - np.array([Gamma_fib(rr) for rr in r]))) / np.max(np.abs(G_init - np.array([Gamma_fib(rr) for rr in r]))):.6f}")
