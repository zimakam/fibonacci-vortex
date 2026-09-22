import math

PHI = (1 + math.sqrt(5))/2
LAM = math.log(PHI)
INV_PHI = 1/PHI
SQRT5 = math.sqrt(5)
FIB = [1,1,2,3,5,8,13,21,34,55,89,144,233,377,610,987,1597,
       2584,4181,6765,10946,17711,28657,46368,75025,121393,
       196418,317811,514229,832040,1346269,2178309,3524578,
       5702887,9227465,14930352,24157817,39088169,63245986,
       102334155,165580141,267914296,433494437,701408733]

# ============================================================
# Method 1: Analytic
# ============================================================
c_analytic = -LAM / (2*(PHI + 2))
print("="*72)
print("METHOD 1: Analytic")
print("="*72)
print()
print(f"  c_inf = -log(phi) / (2*(phi+2))")
print(f"        = -{LAM} / (2 * {PHI+2})")
print(f"        = {c_analytic:.12f}")
print()

# ============================================================
# Method 2: Numerical alpha(K) + linear fit
# ============================================================
def Gamma(r, K):
    return sum(FIB[k]*math.exp(-k*LAM*r) for k in range(K))

def alpha_num(K, h=1e-6):
    g_p = Gamma(1.0 + h, K)
    g_m = Gamma(1.0 - h, K)
    dlog = (math.log(g_p) - math.log(g_m)) / (2*h)
    return -dlog

print("="*72)
print("METHOD 2: Numerical alpha(K), fit diff = c_inf + b/K")
print("="*72)
print()
print(f"{'K':>4} {'alpha_num':>14} {'analytic_part':>14} {'diff':>14} {'diff*K':>14}")
print("-"*72)

data = []
for K in range(5, 40):
    a = alpha_num(K)
    ana = (K-1)*LAM/2 + c_analytic
    diff = a - ana
    data.append((K, a, diff, diff*K))
    if K in (5, 8, 10, 12, 15, 20, 25, 30, 35, 39):
        print(f"{K:>4} {a:>14.8f} {ana:>14.8f} {diff:>14.8f} {diff*K:>14.6f}")

print()
print("If diff -> c_inf, and diff*K -> b:")
print(f"  Last 5 diff: {[f'{d[2]:.6f}' for d in data[-5:]]}")
print(f"  Last 5 diff*K: {[f'{d[3]:.4f}' for d in data[-5:]]}")
print()

# ============================================================
# Method 3: Direct computation of A_inf, B_inf
# ============================================================
print("="*72)
print("METHOD 3: Direct sums A_K, B_K to large K")
print("="*72)
print()
print("A_K = sum_{k=0}^{K-1} (-1)^k * phi^(-2k)")
print("B_K = sum_{k=0}^{K-1} k*(-1)^k * phi^(-2k)")
print()

A_sum = 0.0
B_sum = 0.0
x = INV_PHI**2  # = phi^-2
for K in range(1, 200):
    k = K-1
    A_sum += ((-1)**k) * (x**k)
    B_sum += k * ((-1)**k) * (x**k)

# Closed forms
A_inf = 1 / (1 + x)
B_inf = -x / (1 + x)**2

print(f"  A_inf (numeric sum) = {A_sum:.12f}")
print(f"  A_inf (closed form) = {A_inf:.12f}")
print(f"  A_inf = 1/(1 + phi^-2) = 1/(3-phi) = {1/(3-PHI):.12f}")
print()
print(f"  B_inf (numeric sum) = {B_sum:.12f}")
print(f"  B_inf (closed form) = {B_inf:.12f}")
print(f"  B_inf = -phi^-2 / (1+phi^-2)^2")
print()

# Now compute alpha(K) from analytic expansion
print("="*72)
print("METHOD 3b: alpha(K) from analytic expansion")
print("="*72)
print()
print("Using S1 = (phi*K + phi^-1 * A_K)/sqrt(5)")
print("      S2 = (phi*K(K-1)/2 + phi^-1 * B_K)/sqrt(5)")
print()

def alpha_analytic_full(K):
    A_K = sum(((-1)**k) * (x**k) for k in range(K))
    B_K = sum(k * ((-1)**k) * (x**k) for k in range(K))
    S1 = (PHI*K + INV_PHI*A_K)/SQRT5
    S2 = (PHI*K*(K-1)/2 + INV_PHI*B_K)/SQRT5
    return LAM * S2 / S1

print(f"{'K':>4} {'alpha_num':>14} {'alpha_full_anal':>18} {'diff':>14}")
print("-"*72)
for K in [5, 8, 12, 15, 20, 25, 30, 35, 40]:
    if K > len(FIB):
        continue
    a_num = alpha_num(K)
    a_full = alpha_analytic_full(K)
    print(f"{K:>4} {a_num:>14.8f} {a_full:>18.8f} {a_num - a_full:>14.2e}")

# ============================================================
# Final comparison
# ============================================================
print()
print("="*72)
print("FINAL COMPARISON")
print("="*72)
print()
print(f"  c_inf (analytic formula)   = {c_analytic:.12f}")
print(f"  c_inf (numeric A_K, B_K)   = {-LAM*A_inf/(2*PHI**2):.12f}")
print(f"  Difference                 = {abs(c_analytic + LAM*A_inf/(2*PHI**2)):.2e}")
print()

# Check
alt = -LAM * A_inf / (2 * PHI * PHI)
print(f"  Verify identity:")
print(f"    -log(phi) / (2*(phi+2))             = {c_analytic:.12f}")
print(f"    -log(phi)*A_inf / (2*phi^2)        = {alt:.12f}")
print(f"    A_inf/(1+phi+phi^2) = A_inf/(2+phi) = {A_inf/(2+PHI):.12f}")
print(f"    1/(2+phi) = 1/(phi^2+1)            = {1/(PHI**2+1):.12f}")
print()

# Verify with phi identity
print("  Check: phi^2 = phi + 1, so phi^2 + 1 = phi + 2")
print(f"    phi^2 + 1 = {PHI**2 + 1:.12f}")
print(f"    phi + 2   = {PHI + 2:.12f}")
print(f"    Match: {abs(PHI**2 + 1 - (PHI + 2)):.2e}")
