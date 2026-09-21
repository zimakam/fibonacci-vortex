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
