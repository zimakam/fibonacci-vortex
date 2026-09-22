import numpy as np
import math

PHI = (1+math.sqrt(5))/2
INV_PHI = 1.0/PHI  # = 0.6180339...
LAM_PHI = math.log(PHI)

print("="*72)
print("PROBE #9: Aubry-Andre model at alpha = 1/phi")
print("="*72)
print()
print(f"1/phi = {INV_PHI:.10f}")
print("Hurwitz theorem: 1/phi is the 'most irrational' number")
print()

def aubry_andre_spectrum(N=200, alpha=INV_PHI, lam=1.0, seed=42):
    """Build AA Hamiltonian, return eigenvalues."""
    rng = np.random.default_rng(seed)
    theta = rng.uniform(0, 2*math.pi)
    H = np.zeros((N, N))
    for n in range(N):
        V = 2*lam*math.cos(2*math.pi*alpha*n + theta)
        H[n, n] = V
        if n > 0: H[n, n-1] = 1.0
        if n < N-1: H[n, n+1] = 1.0
    eigs = np.linalg.eigvalsh(H)
    return np.sort(eigs)

# Test at different lambda
print(f"{'lambda':>8} {'N_eig':>6} {'spectrum_range':>18} {'log-density_FFT_peak':>22}")
print("-"*72)

for lam in [0.5, 1.0, 1.5, 2.0]:
    eigs = aubry_andre_spectrum(N=400, alpha=INV_PHI, lam=lam)
    # Density of states -> log-periodic analysis
    # Histogram over energy
    hist, edges = np.histogram(eigs, bins=100, density=True)
    # Find spectrum range
    range_eig = eigs[-1] - eigs[0]
    # Log-periodic: FFT of deviations from smooth
    # Actually simpler: gaps between consecutive eigenvalues
    gaps = np.diff(eigs)
    gaps = gaps[gaps > 1e-6]
    # Distribution of gaps
    log_gaps = np.log(gaps + 1e-12)
    log_gaps -= log_gaps.mean()
    if len(log_gaps) > 10:
        fft = np.abs(np.fft.rfft(log_gaps))
        # Peak excluding DC
        if len(fft) > 2:
            pk = np.argmax(fft[1:]) + 1
            pk_freq = pk
        else:
            pk_freq = 0
    else:
        pk_freq = 0
    print(f"{lam:>8.2f} {len(eigs):>6} {range_eig:>18.4f} {pk_freq:>22}")

print()
print("="*72)
print("Comparison: FibonacciVortex spectrum")
print("="*72)
print()
print("FibonacciVortex has eigenvalues {0, LAM, 2LAM, 3LAM, 4LAM}")
print(f"  step = log(phi) = {LAM_PHI:.6f}")
print()

# Check: does the AA spectral density have phi-related log-periodicity?
# Take lam=1.0 (critical point), compute DOS and its Fourier transform in log E

N = 800
eigs = aubry_andre_spectrum(N=N, alpha=INV_PHI, lam=1.0)
# Normalize energies
e_min, e_max = eigs[0], eigs[-1]
e_norm = (eigs - e_min) / (e_max - e_min)

# Log-periodic analysis of gaps
gaps = np.diff(eigs)
gaps = gaps[gaps > 1e-8]
log_gaps = np.log(gaps)
# Uniform interpolation in log-gap space
log_gaps_sorted = np.sort(log_gaps)
lgs_uniform = np.linspace(log_gaps_sorted.min(), log_gaps_sorted.max(), 256)
# Density
from scipy import stats as sps  # may not have scipy

# Without scipy: simple histogram
hist_lg, edges_lg = np.histogram(log_gaps, bins=64, density=True)
centers = 0.5*(edges_lg[1:] + edges_lg[:-1])
dc = hist_lg - hist_lg.mean()
fft_lg = np.abs(np.fft.rfft(dc))
freqs = np.fft.rfftfreq(len(hist_lg), d=(centers[1]-centers[0]))

print(f"FFT of log-gap density (AA lam=1.0, N={N}):")
print(f"Top-5 peaks:")
top_idx = np.argsort(fft_lg)[::-1][:5]
for i in top_idx:
    print(f"  freq = {freqs[i]:.4f}, amp = {fft_lg[i]:.4f}")

print()
print("Expected phi-related frequency if log-periodic with step log(phi):")
print(f"  f_target = 1/log(phi) = {1/LAM_PHI:.4f}")
print()

# Try also for FibonacciVortex spectrum
print("="*72)
print("Compare with FibonacciVortex gaps")
print("="*72)
fib_gaps = [LAM_PHI]*4  # all gaps = LAM in the SUSY spectrum
print(f"FibonacciVortex eigenvalue gaps: {[f'{g:.4f}' for g in fib_gaps]}")
print("All equal -> no log-periodicity in spectrum itself")
print()
print("But the FUNCTION Gamma_fib(r) DOES have log-periodicity in r.")
print("This is a different aspect.")
