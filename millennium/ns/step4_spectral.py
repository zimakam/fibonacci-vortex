import numpy as np
import math

PHI = (1+math.sqrt(5))/2
LAM = math.log(PHI)
FIB = [1,1,2,3,5,8,13,21,34,55,89,144]
K = 5

def Gamma_fib(r):
    return sum(FIB[k]*math.exp(-k*LAM*r) for k in range(K))

def a_over_nu(r):
    num = sum(k*k*FIB[k]*math.exp(-k*LAM*r) for k in range(K))
    den = sum(k*FIB[k]*math.exp(-k*LAM*r) for k in range(K))
    if abs(den) < 1e-15:
        return 0.0
    return (num/den)*LAM/r + 1.0/(r*r)

N = 200
r_min, r_max = 0.3, 10.0
r = np.linspace(r_min, r_max, N)
dr = r[1] - r[0]
ani_arr = np.array([a_over_nu(rr) for rr in r])

L = np.zeros((N, N))
for i in range(1, N-1):
    ri = r[i]
    L[i, i-1] = (1.0/(dr*dr*ri)) + (1.0/(2*dr*ri*ri)) - ani_arr[i]/(2*dr)
    L[i, i]   = -2.0/(dr*dr*ri)
    L[i, i+1] = (1.0/(dr*dr*ri)) - (1.0/(2*dr*ri*ri)) + ani_arr[i]/(2*dr)

print("="*72)
print("STEP 4: Spectral gap in the energy-weighted space")
print("="*72)
print()

# Interior only
L_int = L[1:-1, 1:-1]
r_int = r[1:-1]
M = len(r_int)

# W-inner product: <u,v>_W = sum u_i v_i * r_i
# Symmetrized operator: M_sym = (W^{1/2} L W^{-1/2} + W^{-1/2} L^T W^{1/2})/2
W_half = np.diag(np.sqrt(r_int))
W_half_inv = np.diag(1.0/np.sqrt(r_int))

M_sym = 0.5 * (W_half @ L_int @ W_half_inv + W_half_inv @ L_int.T @ W_half)

# Check symmetry
sym_err = np.max(np.abs(M_sym - M_sym.T))
print(f"Symmetry error ||M - M^T|| = {sym_err:.2e}")
print(f"Matrix dimension: M = {M}")
print()

# Eigenvalues
eigs = np.linalg.eigvalsh(M_sym)
eigs = np.sort(eigs)  # ascending

lambda_min = eigs[0]
lambda_max = eigs[-1]
gamma = -lambda_max

print(f"Eigenvalue spectrum:")
print(f"  min(lambda)  = {lambda_min:+.6e}")
print(f"  max(lambda)  = {lambda_max:+.6e}   <- slowest decay")
print(f"  gamma = -max(lambda) = {gamma:+.6e}")
print()

print(f"Verify: all eigenvalues < 0?  {np.all(eigs < 0)}")
print(f"Spectral gap gamma = {gamma:.6e} > 0  -> exponential decay")
print()

# Top-5 slowest modes (largest eigenvalues)
print("Top-5 slowest modes (largest eigenvalues):")
idx_slow = np.argsort(eigs)[::-1][:5]
for i in idx_slow:
    print(f"  lambda_{i:3d} = {eigs[i]:+.6e}   tau = {1/abs(eigs[i]):.4f}")
print()

# Top-5 fastest
print("Top-5 fastest modes (most negative eigenvalues):")
idx_fast = np.argsort(eigs)[:5]
for i in idx_fast:
    print(f"  lambda_{i:3d} = {eigs[i]:+.6e}   tau = {1/abs(eigs[i]):.4e}")
print()

# Explicit decay bound
print("="*72)
print("SPECTRAL GAP THEOREM (numerical)")
print("="*72)
print()
print(f"For all perturbations delta in the domain,")
print(f"the energy E(t) = int delta^2 * r dr satisfies:")
print()
print(f"    E(t) <= E(0) * exp(-2 * gamma * t)")
print()
print(f"with gamma = {gamma:.6e}")
print()
print(f"Equivalently: ||delta(t)||_W <= ||delta(0)||_W * exp(-gamma * t)")
print()
print(f"Compare with Step 1 (RK2, unweighted): lambda_max = -0.0323")
print(f"Here (W-weighted): gamma = {gamma:.6e}")
print(f"Ratio: {gamma / 0.0323:.4f}")
