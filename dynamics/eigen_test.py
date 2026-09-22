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

# Build linear operator L such that dG/dt = L G
# L[i] = (G[i+1]-2G[i]+G[i-1])/(dr^2 * r[i]) 
#      - (G[i+1]-G[i-1])/(2 dr * r[i]^2)
#      + ani_arr[i] * (G[i+1]-G[i-1])/(2 dr)
L = np.zeros((N, N))
for i in range(1, N-1):
    L[i, i-1] = 1.0/(dr*dr * r[i]) + 1.0/(2*dr * r[i]*r[i]) - ani_arr[i]/(2*dr)
    L[i, i]   = -2.0/(dr*dr * r[i])
    L[i, i+1] = 1.0/(dr*dr * r[i]) - 1.0/(2*dr * r[i]*r[i]) + ani_arr[i]/(2*dr)
# Boundary: fixed
L[0, 0] = 1.0
L[-1, -1] = 1.0

# Eigenvalues
eigs = np.linalg.eigvals(L)
real_parts = np.real(eigs)
max_real = np.max(real_parts[1:-1])  # skip boundaries
min_real = np.min(real_parts[1:-1])

print("="*60)
print("PROBE #3: Linearized eigenvalues of Burgers operator")
print("="*60)
print()
print(f"Grid: N={N}, r in [{r_min}, {r_max}]")
print()
print(f"Max real part (excluding boundaries): {max_real:+.4e}")
print(f"Min real part: {min_real:+.4e}")
print()
print("Eigenvalues sorted by real part (first 10):")
idx = np.argsort(real_parts)
for i in idx[:10]:
    print(f"  λ_{i:2d} = {eigs[i].real:+.4e} {eigs[i].imag:+.4e}i")

print()
if max_real < 1e-6:
    print("*** All eigenvalues have Re ≤ 0 -> STABLE ***")
elif max_real < 1e-3:
    print("*** Nearly zero real parts -> MARGINAL ***")
else:
    print("*** Positive real parts -> UNSTABLE ***")
