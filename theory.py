import math
from fibonacci_vortex import FibonacciVortex, PHI, INV_PHI
C_INF = -0.0665
LOG_PHI = math.log(PHI)

def alpha_numeric(r, K, h=1e-5):
    fv = FibonacciVortex(k_layers=K)
    gp = fv.gamma_of_r(r + h)
    gm = fv.gamma_of_r(r - h)
    return -(math.log(gp) - math.log(gm)) / (2*h) * r

def alpha_analytic(K):
    return (K - 1) * LOG_PHI / 2 + C_INF

def K_star():
    return 1 + 2*(1 - C_INF) / LOG_PHI

if __name__ == "__main__":
    print("K* = %.4f" % K_star())
    print("K   alpha(r=1)  analytic   diff")
    for K in (5,8,10,12,15,20,23):
        a = alpha_numeric(1.0, K)
        p = alpha_analytic(K)
        print("%2d  %+.6f  %+.6f  %+.6f" % (K, a, p, a-p))
    print()
    print("r   alpha(r) K=5")
    for r in (0.3,0.5,0.8,1.0,1.2,1.5,2.0,3.0):
        print("%4.2f  %+.5f" % (r, alpha_numeric(r, 5)))
    print()
    print("r   Gamma(phi r)/Gamma(r)  1/phi    diff")
    fv = FibonacciVortex(k_layers=5)
    for r in (0.8,1.0,1.2,1.5,2.0):
        ratio = fv.self_similarity(r)
        print("%4.2f  %.6f  %.6f  %+.6f" % (r, ratio, INV_PHI, ratio-INV_PHI))


# ==== Burgers-vortex identity (verified 2026-09-21) ====
def kbar_w(r, K=5, rc=1.0):
    lam = math.log(PHI)/rc
    num = 0.0; den = 0.0
    for k in range(K):
        w = k * FIB[k] * math.exp(-k*lam*r)
        num += k * w
        den += w
    return num/den if den > 0 else 0.0

def a_over_nu(r, K=5, rc=1.0):
    lam = math.log(PHI)/rc
    return kbar_w(r, K, rc)*lam/r + 1.0/r**2

def verify_burgers_identity(r, K=5, rc=1.0, h=1e-5):
    lam = math.log(PHI)/rc
    def dG(r, n):
        return sum((-k*lam)**n * FIB[k]*math.exp(-k*lam*r) for k in range(K))
    lhs = dG(r, 2)/r - dG(r, 1)/r**2
    rhs = -a_over_nu(r, K, rc) * dG(r, 1)
    return lhs - rhs
