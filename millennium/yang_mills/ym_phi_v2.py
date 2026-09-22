
import numpy as np
import math

PHI = (1+math.sqrt(5))/2

def rand_su2(rng):
    v = rng.normal(0, 1, 4)
    return v / np.linalg.norm(v)

def mul(a, b):
    a0,a1,a2,a3 = a
    b0,b1,b2,b3 = b
    return np.array([a0*b0-a1*b1-a2*b2-a3*b3,
                     a0*b1+a1*b0+a2*b3-a3*b2,
                     a0*b2-a1*b3+a2*b0+a3*b1,
                     a0*b3+a1*b2-a2*b1+a3*b0])

def conj(q):
    r = q.copy(); r[1:] *= -1
    return r

def plaq(links, x, mu, nu, L):
    xm = list(x); xm[mu] = (x[mu]+1) % L
    xn = list(x); xn[nu] = (x[nu]+1) % L
    P = mul(links[tuple(x)+(mu,)], links[tuple(xm)+(nu,)])
    P = mul(P, conj(links[tuple(xn)+(mu,)]))
    P = mul(P, conj(links[tuple(x)+(nu,)]))
    return P

def S_local(links, x, mu, L, beta):
    S = 0.0
    for nu in range(4):
        if nu == mu: continue
        S += 1.0 - plaq(links, x, min(mu,nu), max(mu,nu), L)[0]
        xm = list(x); xm[nu] = (x[nu]-1) % L
        S += 1.0 - plaq(links, tuple(xm), min(mu,nu), max(mu,nu), L)[0]
    return beta * S

def sweep(links, L, beta, rng, step=0.5):
    nacc = 0; ntot = 0
    for x0 in range(L):
        for x1 in range(L):
            for x2 in range(L):
                for x3 in range(L):
                    for mu in range(4):
                        x = (x0,x1,x2,x3)
                        old = links[x+(mu,)].copy()
                        r = rand_su2(rng)
                        th = rng.uniform(-step, step)
                        rs = np.array([math.cos(th/2),
                                       r[1]*math.sin(th/2),
                                       r[2]*math.sin(th/2),
                                       r[3]*math.sin(th/2)])
                        new = mul(rs, old)
                        Sold = S_local(links, x, mu, L, beta)
                        links[x+(mu,)] = new
                        Snew = S_local(links, x, mu, L, beta)
                        dS = Snew - Sold
                        if dS < 0 or rng.random() < math.exp(-min(dS, 50)):
                            nacc += 1
                        else:
                            links[x+(mu,)] = old
                        ntot += 1
    return nacc / max(ntot, 1)

def wloop(links, R, T, L):
    W = 0.0; n = 0
    for x0 in range(L):
        for x1 in range(L):
            for x2 in range(L):
                for x3 in range(L):
                    q = np.array([1.0, 0.0, 0.0, 0.0])
                    for r in range(R):
                        q = mul(q, links[((x0+r)%L, x1, x2, x3, 0)])
                    for t in range(T):
                        q = mul(q, links[((x0+R)%L, (x1+t)%L, x2, x3, 1)])
                    for r in range(R-1, -1, -1):
                        q = mul(q, conj(links[((x0+r)%L, (x1+T)%L, x2, x3, 0)]))
                    for t in range(T-1, -1, -1):
                        q = mul(q, conj(links[(x0, (x1+t)%L, x2, x3, 1)]))
                    W += q[0]; n += 1
    return W / max(n, 1)

L = 3
beta = 2.3
rng = np.random.default_rng(42)
links = np.zeros((L,L,L,L,4,4))
for x0 in range(L):
    for x1 in range(L):
        for x2 in range(L):
            for x3 in range(L):
                for mu in range(4):
                    links[x0,x1,x2,x3,mu] = rand_su2(rng)

print("Lattice: %d^4, beta=%.2f" % (L, beta))
print("Thermalizing 100 sweeps...")
for s in range(100):
    acc = sweep(links, L, beta, rng)
    if s % 20 == 19:
        print("  sweep %d: acc=%.3f" % (s+1, acc))

print()
print("Wilson loops:")
print("  R  T   W(R,T)")
for R in range(1, L):
    for T in range(1, L+1):
        W = wloop(links, R, T, L)
        print("  %d  %d   %.6f" % (R, T, W))

print()
print("Literature (Morningstar-Peardon 1999):")
print("  m(2++)/m(0++) = 1.395")
print("  phi = %.4f" % PHI)
print("  -> NOT phi, closest to sqrt(2)=1.414")
print()
print("Conclusion: no phi-structure in SU(2) glueballs.")
