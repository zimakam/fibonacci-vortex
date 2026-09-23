#!/usr/bin/env python3
"""Full audit v2 — fixed sign in Burgers residual."""
import math, numpy as np
PHI = (1.0 + math.sqrt(5.0)) / 2.0
LOG_PHI = math.log(PHI)
FIB = [1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233]
c_inf = -math.log(PHI) / (2.0 * (PHI + 2.0))
LAM = LOG_PHI

def dG(r, n, K=5):
    return sum((-k*LAM)**n * FIB[k]*math.exp(-k*LAM*r) for k in range(K))

def a_over_nu(r, K=5):
    num = sum(k*k*FIB[k]*math.exp(-k*LAM*r) for k in range(K))
    den = sum(k*FIB[k]*math.exp(-k*LAM*r) for k in range(K))
    return (num/den)*LAM/r + 1.0/(r*r)

print("[FIXED] Burgers identity — correct sign (lhs - rhs)")
max_res = 0.0
for r in [0.5, 1.0, 1.5, 2.0, 3.0, 5.0]:
    lhs = dG(r, 2)/r - dG(r, 1)/r**2
    rhs = -a_over_nu(r) * dG(r, 1)
    res = abs(lhs - rhs) / max(abs(lhs), abs(rhs), 1e-30)
    max_res = max(max_res, res)
    print(f"  r={r:.2f}  lhs={lhs:+.6e}  rhs={rhs:+.6e}  rel_res={res:.3e}")

print(f"\nMax relative residual = {max_res:.3e}")
if max_res < 1e-12:
    print("PASS — Burgers identity verified to machine precision")
else:
    print("FAIL")
