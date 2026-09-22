import numpy as np
import math

PHI = (1+math.sqrt(5))/2
LAM = math.log(PHI)
FIB = [1,1,2,3,5,8,13,21,34,55,89,144]
K = 5

def a_over_nu(r):
    num = sum(k*k*FIB[k]*math.exp(-k*LAM*r) for k in range(K))
    den = sum(k*FIB[k]*math.exp(-k*LAM*r) for k in range(K))
    if abs(den) < 1e-15: return 0.0
    return (num/den)*LAM/r + 1.0/(r*r)

N = 80
r_min, r_max = 0.5, 6.0
r = np.linspace(r_min, r_max, N)
dr = r[1] - r[0]
ani_arr = np.array([a_over_nu(rr) for rr in r])

# Interior only: M = N-2 points
M = N - 2
L = np.zeros((M, M))
for i in range(M):
    ii = i + 1  # interior index
    if i > 0:
        L[i, i-1] = 1.0/(dr*dr * r[ii]) + 1.0/(2*dr * r[ii]*r[ii]) - ani_arr[ii]/(2*dr)
    L[i, i]   = -2.0/(dr*dr * r[ii])
    if i + 1 < M:
        L[i, i+1] = 1.0/(dr*dr * r[ii]) - 1.0/(2*dr * r[ii]*r[ii]) + ani_arr[ii]/(2*dr)

eigs = np.linalg.eigvals(L)
real_parts = np.real(eigs)
max_real = np.max(real_parts)
min_real = np.min(real_parts)

print("="*60)
print("PROBE #3 (fixed): Interior eigenvalues only")
print("="*60)
print()
print(f"Interior points: M={M}")
print(f"Max real part: {max_real:+.6e}")
print(f"Min real part: {min_real:+.6e}")
print()
print("Top-10 eigenvalues by real part (largest first):")
idx = np.argsort(real_parts)[::-1]
for i in idx[:10]:
    print(f"  λ_{i:2d} = {eigs[i].real:+.6e} {eigs[i].imag:+.6e}i")

print()
if max_real < -1e-6:
    print("*** All Re(λ) < 0 -> STABLE ***")
    print("*** Perturbations decay exponentially ***")
elif abs(max_real) < 1e-6:
    print("*** All Re(λ) ≈ 0 -> MARGINAL ***")
    print("*** Perturbations neither grow nor decay ***")
else:
    print(f"*** Some Re(λ) > 0 -> UNSTABLE ***")
    print(f"*** Max growth rate: {max_real:+.4e} ***")
