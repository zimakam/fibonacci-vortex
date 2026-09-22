import numpy as np
import math

PHI = (1+math.sqrt(5))/2
LAM = math.log(PHI)

# ============================================================
# Lattice SU(2) Yang-Mills (small lattice)
# ============================================================

def random_su2(rng):
    """Random SU(2) element as quaternion (a0, a1, a2, a3)."""
    v = rng.normal(0, 1, 4)
    v /= np.linalg.norm(v)
    return v

def su2_mul(a, b):
    """Quaternion multiplication for SU(2)."""
    a0, a1, a2, a3 = a
    b0, b1, b2, b3 = b
    return np.array([
        a0*b0 - a1*b1 - a2*b2 - a3*b3,
        a0*b1 + a1*b0 + a2*b3 - a3*b2,
        a0*b2 - a1*b3 + a2*b0 + a3*b1,
        a0*b3 + a1*b2 - a2*b1 + a3*b0,
    ])

def su2_conj(q):
    out = q.copy()
    out[1:] *= -1
    return out

def plaquette(links, x, mu, nu, L):
    """Plaquette at site x in mu-nu plane."""
    x_mu = list(x); x_mu[mu] = (x[mu]+1) % L
    x_nu = list(x); x_nu[nu] = (x[nu]+1) % L
    U1 = links[tuple(x) + (mu,)]
    U2 = links[tuple(x_mu) + (nu,)]
    U3 = links[tuple(x_nu) + (mu,)]
    U4 = links[tuple(x) + (nu,)]
    P = su2_mul(U1, U2)
    P = su2_mul(P, su2_conj(U3))
    P = su2_mul(P, su2_conj(U4))
    return P

def total_action(links, L, beta):
    S = 0.0
    for x0 in range(L):
        for x1 in range(L):
            for x2 in range(L):
                for x3 in range(L):
                    for mu in range(4):
                        for nu in range(mu+1, 4):
                            P = plaquette(links, (x0,x1,x2,x3), mu, nu, L)
                            S += 1.0 - P[0]  # real part
    return beta * S

def local_action(links, x, mu, L, beta):
    """Sum of plaquettes containing link (x, mu)."""
    S = 0.0
    for nu in range(4):
        if nu == mu: continue
        # Plaquette at x
        P = plaquette(links, x, min(mu,nu), max(mu,nu), L)
        S += 1.0 - P[0]
        # Plaquette at x-nu
        x_minus_nu = list(x); x_minus_nu[nu] = (x[nu]-1) % L
        P2 = plaquette(links, tuple(x_minus_nu), min(mu,nu), max(mu,nu), L)
        S += 1.0 - P2[0]
    return beta * S

def metropolis_sweep(links, L, beta, rng, step=0.3):
    """One Metropolis sweep."""
    n_acc = 0
    n_tot = 0
    for x0 in range(L):
        for x1 in range(L):
            for x2 in range(L):
                for x3 in range(L):
                    for mu in range(4):
                        x = (x0, x1, x2, x3)
                        old = links[x + (mu,)]
                        # Random SU(2) proposal
                        r = random_su2(rng)
                        # Rotate by small angle
                        theta = rng.uniform(-step, step)
                        r_small = np.array([math.cos(theta/2),
                                            r[1]*math.sin(theta/2),
                                            r[2]*math.sin(theta/2),
                                            r[3]*math.sin(theta/2)])
                        new = su2_mul(r_small, old)
                        
                        # Delta action
                        links[x + (mu,)] = old
                        S_old = local_action(links, x, mu, L, beta)
                        links[x + (mu,)] = new
                        S_new = local_action(links, x, mu, L, beta)
                        dS = S_new - S_old
                        
                        if dS < 0 or rng.random() < math.exp(-dS):
                            n_acc += 1
                        else:
                            links[x + (mu,)] = old
                        n_tot += 1
    return n_acc / max(n_tot, 1)

def wilson_loop(links, R, T, L):
    """Average Wilson loop of size R x T."""
    W_sum = 0.0
    n = 0
    for x0 in range(L):
        for x1 in range(L):
            for x2 in range(0, min(2,L)):
                for x3 in range(0, min(2,L)):
                    q = np.array([1.0, 0, 0, 0])
                    # Path: (x, 0) R times, (x+R, 1) T times
                    # (x+R, 0) R times reversed, (x, 1) T times reversed
                    for r in range(R):
                        q = su2_mul(q, links[((x0+r)%L, x1, x2, x3, 0)])
                    for t in range(T):
                        q = su2_mul(q, links[((x0+R)%L, (x1+t)%L, x2, x3, 1)])
                    for r in range(R-1, -1, -1):
                        q = su2_mul(q, su2_conj(links[((x0+r)%L, (x1+T)%L, x2, x3, 0)]))
                    for t in range(T-1, -1, -1):
                        q = su2_mul(q, su2_conj(links[(x0, (x1+t)%L, x2, x3, 1)]))
                    W_sum += q[0]
                    n += 1
    return W_sum / max(n, 1)

# ============================================================
# Main computation
# ============================================================
print("="*72)
print("YANG-MILLS: phi-structure in SU(2) lattice spectrum")
print("="*72)
print()

L = 3   # small lattice for speed
beta = 2.3  # standard critical region
rng = np.random.default_rng(42)

# Initialize links
links = np.zeros((L, L, L, L, 4, 4))
for x0 in range(L):
    for x1 in range(L):
        for x2 in range(L):
            for x3 in range(L):
                for mu in range(4):
                    links[x0,x1,x2,x3,mu] = random_su2(rng)

print(f"Lattice: {L}^4, beta={beta}")
print(f"Total links: {L**4 * 4} = {L**4 * 4}")
print()

# Thermalize
print("Thermalizing (20 sweeps)...")
for sweep in range(20):
    acc = metropolis_sweep(links, L, beta, rng)
    if sweep % 5 == 4:
        S = total_action(links, L, beta)
        print(f"  sweep {sweep+1}: accept={acc:.3f}, S={S:.4f}")
print()

# Wilson loops
print("Measuring Wilson loops:")
print(f"{'R':>3} {'T':>3} {'<W(R,T)>':>12} {'ln(W)':>12}")
print("-"*34)
Ws = {}
for R in range(1, min(3, L)):
    for T in range(1, min(4, L)):
        W = wilson_loop(links, R, T, L)
        Ws[(R,T)] = W
        lnW = math.log(W) if W > 0 else -50
        print(f"{R:>3} {T:>3} {W:>12.6f} {lnW:>12.6f}")

# Estimate mass gap from W ~ exp(-m*T) for fixed R
print()
print("Mass gap estimate (slope of ln W vs T, fixed R):")
masses = []
for R in range(1, min(3, L)):
    Ts = [T for T in range(1, min(4, L)) if (R,T) in Ws]
    lnWs = [math.log(Ws[(R,T)]) if Ws[(R,T)] > 0 else -50 for T in Ts]
    if len(Ts) >= 2:
        A = np.vstack([Ts, np.ones_like(Ts)]).T
        try:
            slope, intercept = np.linalg.lstsq(A, lnWs, rcond=None)[0]
            m = -slope
            masses.append(m)
            print(f"  R={R}: m = {m:.6f}")
        except Exception:
            pass

# Phi test: ratios of masses
if len(masses) >= 2:
    print()
    print("Ratios of masses (test for phi):")
    for i in range(len(masses)-1):
        ratio = masses[i+1] / masses[i] if masses[i] > 0 else 0
        print(f"  m_{i+1}/m_{i} = {ratio:.6f}   phi = {PHI:.6f}")

print()
print("="*72)
print("INTERPRETATION")
print("="*72)
print()
print("Standard result (literature): SU(2) glueball masses")
print("  0++ mass / sqrt(sigma) ~ 3.5, no phi-structure.")
print()
print("Our test (very small L=3, beta=2.3): masses are dominated by")
print("lattice artifacts. Not reliable for global structure.")
print()
print("Expected result: no phi-structure in SU(2) spectrum.")
print("Reason: QCD scale is set by dimensional transmutation,")
print("not by Fibonacci/golden ratio.")
