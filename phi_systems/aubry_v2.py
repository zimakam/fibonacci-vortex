import numpy as np
import math

PHI = (1+math.sqrt(5))/2
INV_PHI = 1/PHI
LAM_PHI = math.log(PHI)

def aubry_spectrum(N=2000, alpha=INV_PHI, lam=1.0, theta=0.0):
    H = np.zeros((N, N))
    for n in range(N):
        V = 2*lam*math.cos(2*math.pi*alpha*n + theta)
        H[n, n] = V
        if n > 0: H[n, n-1] = 1.0
        if n < N-1: H[n, n+1] = 1.0
    return np.sort(np.linalg.eigvalsh(H))

def log_periodic_dos(eigs, nbins=200):
    """Density of states, then Fourier over log|E-E_F| near center."""
    # Energy gaps
    gaps = np.diff(eigs)
    gaps = gaps[gaps > 1e-10]
    log_gaps = np.log(gaps)
    # Histogram
    hist, edges = np.histogram(log_gaps, bins=nbins, density=True)
    centers = 0.5*(edges[1:] + edges[:-1])
    dlog = centers[1] - centers[0]
    dc = hist - hist.mean()
    fft = np.abs(np.fft.rfft(dc))
    freqs = np.fft.rfftfreq(len(hist), d=dlog)
    return freqs, fft, centers, hist

print("="*72)
print("PROBE #9 v2: Aubry-Andre with multiple phases + larger N")
print("="*72)
print()

# Average over phases
N = 1500
phases = [0, math.pi/4, math.pi/2, 3*math.pi/4, math.pi, 5*math.pi/4, 3*math.pi/2, 7*math.pi/4]
# Or random
rng = np.random.default_rng(42)
phases = rng.uniform(0, 2*math.pi, 20)

print(f"N={N}, averaging over {len(phases)} random phases")
print()

# Average FFT over phases
def averaged_fft(alpha, lam, phases, N=1500, nbins=200):
    fft_sum = None
    freqs = None
    for theta in phases:
        eigs = aubry_spectrum(N=N, alpha=alpha, lam=lam, theta=theta)
        f, fft, _, _ = log_periodic_dos(eigs, nbins=nbins)
        if fft_sum is None:
            fft_sum = fft
            freqs = f
        else:
            fft_sum = fft_sum + fft
    return freqs, fft_sum / len(phases)

# Test critical point
print("Computing at critical point (lam=1.0)...")
freqs, fft_avg = averaged_fft(INV_PHI, 1.0, phases, N=N)

# Top peaks
idx = np.argsort(fft_avg[1:])[::-1][:8] + 1
print(f"Top-8 peaks (α=1/φ, λ=1.0):")
for i in idx:
    print(f"  freq = {freqs[i]:.4f}, amp = {fft_avg[i]:.4f}")

print()
print(f"f_target = 1/log(phi) = {1/LAM_PHI:.4f}")
print(f"2*f_target = {2/LAM_PHI:.4f}")
print(f"3*f_target = {3/LAM_PHI:.4f}")

# Compare with sqrt(2)-1
print()
print("Computing for α = sqrt(2)-1...")
freqs2, fft_avg2 = averaged_fft(math.sqrt(2)-1, 1.0, phases, N=N)
idx2 = np.argsort(fft_avg2[1:])[::-1][:8] + 1
print(f"Top-8 peaks (α=√2−1, λ=1.0):")
for i in idx2:
    print(f"  freq = {freqs2[i]:.4f}, amp = {fft_avg2[i]:.4f}")

# Direct comparison of spectra
print()
print("="*72)
print("Direct spectrum comparison")
print("="*72)
eigs_phi = aubry_spectrum(N=800, alpha=INV_PHI, lam=1.0, theta=0)
eigs_sq2 = aubry_spectrum(N=800, alpha=math.sqrt(2)-1, lam=1.0, theta=0)
print(f"φ spectrum: min={eigs_phi[0]:.4f}, max={eigs_phi[-1]:.4f}, range={eigs_phi[-1]-eigs_phi[0]:.4f}")
print(f"√2-1 spectrum: min={eigs_sq2[0]:.4f}, max={eigs_sq2[-1]:.4f}, range={eigs_sq2[-1]-eigs_sq2[0]:.4f}")
print()

# Check for exact degeneracies/self-similar structure
gaps_phi = np.diff(eigs_phi)
gaps_phi = gaps_phi[gaps_phi > 1e-8]
print(f"φ gaps: n={len(gaps_phi)}, mean={gaps_phi.mean():.4f}, std={gaps_phi.std():.4f}")
print(f"      : min={gaps_phi.min():.6f}, max={gaps_phi.max():.6f}")
# Check for Fibonacci-like ratios in gaps
ratios = gaps_phi[1:] / gaps_phi[:-1]
# Find φ-closest ratios
phi_close = np.sum(np.abs(ratios - PHI) < 0.01)
print(f"      : ratios close to φ (within 1%): {phi_close} of {len(ratios)}")
