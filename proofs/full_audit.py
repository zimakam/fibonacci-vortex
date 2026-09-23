#!/usr/bin/env python3
"""
Full verification of all claims — FibonacciVortex + Merkaba framework.
Author: Ziyavutdinov M.K. (Zimaka), ORCID 0009-0005-9212-9921
Date: 2026-09-23
Run: python3 full_audit.py
"""
import math
import numpy as np

PHI = (1.0 + math.sqrt(5.0)) / 2.0
INV_PHI = PHI - 1.0
LOG_PHI = math.log(PHI)
FIB = [1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610]

passed = 0
failed = 0
def check(name, cond, detail=""):
    global passed, failed
    if cond:
        print(f"  PASS  {name}  {detail}")
        passed += 1
    else:
        print(f"  FAIL  {name}  {detail}")
        failed += 1

print("=" * 72)
print("FULL AUDIT — verification of all claims")
print("=" * 72)

# ────────────────────────────────────────────────────────────
print("\n[1] c_inf exact formula")
c_inf_exact = -math.log(PHI) / (2.0 * (PHI + 2.0))
c_inf_numeric = -0.066501838644
check("c_inf = -log(phi)/(2(phi+2))",
      abs(c_inf_exact - c_inf_numeric) < 1e-12,
      f"value = {c_inf_exact:.15f}")

# ────────────────────────────────────────────────────────────
print("\n[2] K* = 5.43 and K=5 selection")
K_star = 1 + 2*(1 - c_inf_exact) / LOG_PHI
check("K* = 5.43 (within 0.01)", abs(K_star - 5.43) < 0.01,
      f"K* = {K_star:.6f}")
check("K* closer to 5 than 6", abs(K_star - 5) < abs(K_star - 6),
      f"|K*-5|={abs(K_star-5):.4f}  |K*-6|={abs(K_star-6):.4f}")

# ────────────────────────────────────────────────────────────
print("\n[3] S1(0) = sum k*F_k = 34 = F_9")
K = 5
S1_at_0 = sum(k * FIB[k] for k in range(K))
check("sum k*F_k for k=0..4 equals 34",
      S1_at_0 == 34, f"= {S1_at_0}")
check("34 = F_9 (FIB[8]=34)", FIB[8] == 34)

# ────────────────────────────────────────────────────────────
print("\n[4] Identity 1+2+6+12 = 21 = F_8")
deg_sum = 1 + 2 + 6 + 12
check("1+2+6+12 = 21", deg_sum == 21)
check("21 = F_8 (FIB[7]=21)", FIB[7] == 21)

# ────────────────────────────────────────────────────────────
print("\n[5] FIB_SPEEDS = (21, 13, 8, 5)")
FIB_SPEEDS = (21.0, 13.0, 8.0, 5.0)
check("FIB_SPEEDS == (F_8, F_7, F_6, F_5)",
      FIB_SPEEDS == (FIB[7], FIB[6], FIB[5], FIB[4]),
      f"= {FIB_SPEEDS}")
ratios = [FIB_SPEEDS[i] / FIB_SPEEDS[i+1] for i in range(3)]
print(f"     ratios: {[f'{r:.4f}' for r in ratios]}  (phi={PHI:.4f})")
check("all ratios close to phi (within 0.03)",
      all(abs(r - PHI) < 0.03 for r in ratios))

# ────────────────────────────────────────────────────────────
print("\n[6] R_ext / R_int = phi^3")
N_PAIRS = 4
R_ext_over_R_int = PHI ** (N_PAIRS - 1)
check("phi^3 = 4.236", abs(R_ext_over_R_int - 4.2360679) < 1e-6,
      f"phi^3 = {R_ext_over_R_int:.8f}")

# ────────────────────────────────────────────────────────────
print("\n[7] lambda_Merkaba = omega^2 * sqrt(5) / c^2")
omega_0 = 1.0
c = 10.0
ou = +omega_0 * PHI
ol = -omega_0 * INV_PHI
lam_from_omegas = (ou**2 - ol**2) / (c**2)
lam_formula = (omega_0**2 * math.sqrt(5.0)) / (c**2)
check("(ou^2 - ol^2)/c^2 == omega^2 * sqrt(5) / c^2",
      abs(lam_from_omegas - lam_formula) < 1e-15,
      f"diff = {abs(lam_from_omegas - lam_formula):.3e}")

# golden ratio identity phi^2 - 1/phi^2 = sqrt(5)
phi_sq_minus = PHI**2 - INV_PHI**2
check("phi^2 - 1/phi^2 = sqrt(5)",
      abs(phi_sq_minus - math.sqrt(5.0)) < 1e-14,
      f"diff = {abs(phi_sq_minus - math.sqrt(5.0)):.3e}")

# ────────────────────────────────────────────────────────────
print("\n[8] alpha(K) asymptotic formula")
def alpha_analytic(K):
    return (K - 1) * LOG_PHI / 2.0 + c_inf_exact

alpha_5 = alpha_analytic(5)
alpha_8 = alpha_analytic(8)
print(f"     alpha(5) = {alpha_5:+.6f}")
print(f"     alpha(8) = {alpha_8:+.6f}")
check("alpha(5) crosses 1 (alpha(5) slightly less than 1)",
      alpha_5 < 1.0 and alpha_5 > 0.8,
      f"alpha(5) = {alpha_5:.6f}")
check("alpha(8) > 1", alpha_8 > 1.0)

# ────────────────────────────────────────────────────────────
print("\n[9] Burgers identity for Gamma(r) = sum F_k phi^(-k LAM r)")
LAM = LOG_PHI

def Gamma_r(r, K=5):
    return sum(FIB[k] * math.exp(-k * LAM * r) for k in range(K))

def dG(r, n, K=5):
    return sum((-k*LAM)**n * FIB[k] * math.exp(-k*LAM*r) for k in range(K))

def a_over_nu(r, K=5):
    num = sum(k*k*FIB[k]*math.exp(-k*LAM*r) for k in range(K))
    den = sum(k*FIB[k]*math.exp(-k*LAM*r) for k in range(K))
    if abs(den) < 1e-15: return 0.0
    return (num/den)*LAM/r + 1.0/(r*r)

max_residual = 0.0
for r in [0.5, 1.0, 1.5, 2.0, 3.0, 5.0]:
    lhs = dG(r, 2)/r - dG(r, 1)/r**2
    rhs = -a_over_nu(r) * dG(r, 1)
    residual = abs(lhs + rhs) / max(abs(lhs), abs(rhs), 1e-30)
    max_residual = max(max_residual, residual)
check("Burgers identity residual < 1e-12",
      max_residual < 1e-12,
      f"max residual = {max_residual:.3e}")

# ────────────────────────────────────────────────────────────
print("\n[10] Roots of S1 polynomial (K=5)")
coef = [k * FIB[k] for k in range(K)]
p = coef[::-1]
roots = np.roots(p)
print(f"     S1 coeffs (k*F_k, k=0..4): {coef}")
print(f"     roots of S1 (in x = e^(-LAM r)):")
for rt in roots:
    print(f"       x = {rt.real:+.6f} {rt.imag:+.6f}i   |x|={abs(rt):.6f}  arg={np.degrees(np.angle(rt)):+.3f}°")
# trivial root at 0
nontrivial = [rt for rt in roots if abs(rt) > 1e-10]
check("S1 has 3 nontrivial roots", len(nontrivial) == 3,
      f"found {len(nontrivial)}")
# conjugate pair + real root
has_real = any(abs(rt.imag) < 1e-10 for rt in nontrivial)
has_complex = any(abs(rt.imag) > 1e-10 for rt in nontrivial)
check("1 real + 1 conjugate pair", has_real and has_complex)

# ────────────────────────────────────────────────────────────
print("\n[11] Mirror theorem Gamma(r, pi-theta) = Gamma_alt(r, theta)")
def Gamma_theta(r, theta, K=5, alt=False):
    s = 0.0
    for k in range(K):
        sign = (-1)**k if alt else 1
        s += sign * FIB[k] * math.exp(-k*LAM*r) * math.cos(k*theta)
    return s

max_mirror_err = 0.0
for r in [0.5, 1.0, 2.0]:
    for theta in [0.3, 1.0, 1.5, 2.0]:
        lhs = Gamma_theta(r, math.pi - theta, alt=False)
        rhs = Gamma_theta(r, theta, alt=True)
        max_mirror_err = max(max_mirror_err, abs(lhs - rhs))
check("mirror theorem residual < 1e-14",
      max_mirror_err < 1e-14,
      f"max err = {max_mirror_err:.3e}")

# ────────────────────────────────────────────────────────────
print("\n[12] Schottky anomaly at beta* = 2.04")
def Z(beta, K=5):
    return sum(k * FIB[k] * math.exp(-k*beta) for k in range(K) if k > 0)

def U(beta, K=5):
    num = sum(k*k * FIB[k] * math.exp(-k*beta) for k in range(K) if k > 0)
    den = Z(beta, K)
    return num / den if den > 0 else 0.0

def C(beta, K=5, h=1e-5):
    dU = (U(beta+h) - U(beta-h)) / (2*h)
    return -beta*beta * dU

# scan for maximum of C
betas = np.linspace(0.01, 5.0, 2000)
Cs = np.array([C(b) for b in betas])
i_max = np.argmax(Cs)
beta_star = betas[i_max]
C_max = Cs[i_max]
print(f"     beta* = {beta_star:.4f}   C_max = {C_max:.4f}")
check("beta* near 2.04 (within 0.3)",
      abs(beta_star - 2.04) < 0.3,
      f"beta* = {beta_star:.4f}")
check("C_max positive and bounded",
      0 < C_max < 5.0,
      f"C_max = {C_max:.4f}")

# ────────────────────────────────────────────────────────────
print("\n[13] DNA mirror structure j = n-1-i")
n_genes = 21
cross_weights = [PHI**(-min(i, n_genes-1-i)/n_genes) for i in range(n_genes)]
check("21 genes = F_8", n_genes == 21)
check("weights symmetric: w[i] == w[j]",
      all(abs(cross_weights[i] - cross_weights[n_genes-1-i]) < 1e-15
          for i in range(n_genes)))
check("central weight w[10] = phi^(-10/21)",
      abs(cross_weights[10] - PHI**(-10.0/21.0)) < 1e-15,
      f"w[10] = {cross_weights[10]:.6f}")

# ────────────────────────────────────────────────────────────
print()
print("=" * 72)
print(f"RESULT: {passed} PASSED, {failed} FAILED")
print("=" * 72)
