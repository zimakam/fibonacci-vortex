
import numpy as np
import math

PHI = (1+math.sqrt(5))/2

def mul(a, b):
    a0,a1,a2,a3 = a[...,0], a[...,1], a[...,2], a[...,3]
    b0,b1,b2,b3 = b[...,0], b[...,1], b[...,2], b[...,3]
    return np.stack([
        a0*b0 - a1*b1 - a2*b2 - a3*b3,
        a0*b1 + a1*b0 + a2*b3 - a3*b2,
        a0*b2 - a1*b3 + a2*b0 + a3*b1,
        a0*b3 + a1*b2 - a2*b1 + a3*b0], axis=-1)

def conj(q):
    r = q.copy()
    r[..., 1:] *= -1
    return r

def plaquette(links, x, mu, nu, L):
    xm = list(x); xm[mu] = (x[mu]+1) % L
    xn = list(x); xn[nu] = (x[nu]+1) % L
    U1 = links[tuple(x)+(mu,)]
    U2 = links[tuple(xm)+(nu,)]
    U3 = links[tuple(xn)+(mu,)]
    U4 = links[tuple(x)+(nu,)]
    P = mul(U1, U2)
    P = mul(P, conj(U3))
    P = mul(P, conj(U4))
    return P

def staple_sum(links, mu):
    """K_staple such that (1/2) Tr[U_mu * K_staple] = sum plaquette traces containing U_mu."""
    L = links.shape[0]
    K = np.zeros_like(links[..., mu, :])
    for nu in range(4):
        if nu == mu:
            continue
        # K_f_rev = U_nu(x) * U_mu(x+nu) * U_nu^dag(x+mu)
        U_nu_x = links[..., nu, :]
        U_mu_xnu = np.roll(links[..., mu, :], -1, axis=nu)
        U_nu_xmu = np.roll(links[..., nu, :], -1, axis=mu)
        K_f = mul(U_nu_x, U_mu_xnu)
        K_f = mul(K_f, conj(U_nu_xmu))
        K = K + K_f
        # K_b = U_nu^dag(x-nu) * U_mu(x-nu) * U_nu(x-nu+mu)
        U_nu_xmnu = np.roll(links[..., nu, :], 1, axis=nu)
        U_mu_xmnu = np.roll(links[..., mu, :], 1, axis=nu)
        U_nu_xmnu_xmu = np.roll(U_nu_xmnu, -1, axis=mu)
        K_b = mul(conj(U_nu_xmnu), U_mu_xmnu)
        K_b = mul(K_b, U_nu_xmnu_xmu)
        K = K + K_b
    return K

def verify_staple(links, beta):
    """Check S_direct = S_staple/2."""
    L = links.shape[0]
    S_direct = 0.0
    for x0 in range(L):
        for x1 in range(L):
            for x2 in range(L):
                for x3 in range(L):
                    for mu in range(4):
                        for nu in range(mu+1, 4):
                            P = plaquette(links, (x0,x1,x2,x3), mu, nu, L)
                            S_direct += 1.0 - P[0]
    S_direct *= beta

    S_staple = 0.0
    for mu in range(4):
        K = staple_sum(links, mu)
        U = links[..., mu, :]
        # dot product U·K (4-vector)
        dot = np.sum(U * K, axis=-1)
        S_staple += np.sum(6*beta - beta*dot)
    return S_direct, S_staple/2

def sweep(links, beta, rng, step=1.5):
    L = links.shape[0]
    nacc = 0; ntot = 0
    for mu in range(4):
        K = staple_sum(links, mu)
        r = rng.normal(0, 1, (L,L,L,L,4))
        r /= np.linalg.norm(r, axis=-1, keepdims=True)
        th = rng.uniform(-step, step, (L,L,L,L))
        rs = np.stack([
            np.cos(th/2),
            r[...,1]*np.sin(th/2),
            r[...,2]*np.sin(th/2),
            r[...,3]*np.sin(th/2)], axis=-1)
        old = links[..., mu, :].copy()
        new = mul(rs, old)
        diff = new - old
        # CORRECT: dS = -beta * (U' - U) . K_staple
        dS = -beta * np.sum(diff * K, axis=-1)
        acc_prob = np.where(dS < 0, 1.0, np.exp(-np.minimum(dS, 50)))
        accept = rng.random((L,L,L,L)) < acc_prob
        links[..., mu, :] = np.where(accept[..., None], new, old)
        nacc += np.sum(accept)
        ntot += L**4
    return nacc / max(ntot, 1)

def wloop(links, R, T):
    L = links.shape[0]
    q = np.zeros((L,L,L,L,4)); q[..., 0] = 1.0
    for r in range(R):
        q = mul(q, np.roll(links[..., 0, :], -r, axis=0))
    for t in range(T):
        s = np.roll(links[..., 1, :], -t, axis=1)
        s = np.roll(s, -R, axis=0)
        q = mul(q, s)
    for r in range(R-1, -1, -1):
        s = np.roll(links[..., 0, :], -r, axis=0)
        s = np.roll(s, -T, axis=1)
        q = mul(q, conj(s))
    for t in range(T-1, -1, -1):
        q = mul(q, conj(np.roll(links[..., 1, :], -t, axis=1)))
    return float(np.mean(q[..., 0]))

L = 4
beta = 2.3
rng = np.random.default_rng(42)
links = rng.normal(0, 1, (L,L,L,L,4,4))
links /= np.maximum(np.linalg.norm(links, axis=-1, keepdims=True), 1e-12)

# Verify staple formula
S_d, S_s = verify_staple(links, beta)
print("Verify staple: S_direct=%.4f  S_staple/2=%.4f  diff=%.2e" % (S_d, S_s, abs(S_d - S_s)))
print()

print("Lattice: %d^4, beta=%.2f" % (L, beta))
print("Thermalizing 200 sweeps (step=1.5)...")
for s in range(200):
    acc = sweep(links, beta, rng, step=1.5)
    if s % 50 == 49:
        print("  sweep %d: acc=%.3f" % (s+1, acc))

print()
print("Wilson loops:")
print("  R  T   W")
for R in range(1, L):
    for T in range(1, L+1):
        W = wloop(links, R, T)
        print("  %d  %d   %.6f" % (R, T, W))

print()
print("phi = %.4f" % PHI)
print("sqrt(2) = %.4f" % math.sqrt(2))
print("Literature: m(2++)/m(0++) = 1.395 (closest to sqrt(2)))")
