import numpy as np
import math

PHI = (1+math.sqrt(5))/2
INV_PHI = 1.0/PHI
LAM_PHI = math.log(PHI)

print("="*72)
print("PROBE #9 (fixed): Aubry-Andre log-periodic structure")
print("="*72)
print()
print(f"1/phi = {INV_PHI:.10f}")
print(f"f_target = 1/log(phi) = {1/LAM_PHI:.6f}")
print()

def aubry_spectrum(N=800, alpha=INV_PHI, lam=1.0, seed=42):
    rng = np.random.default_rng(seed)
    theta = rng.uniform(0, 2*math.pi)
    H = np.zeros((N, N))
    for n in range(N):
        V = 2*lam*math.cos(2*math.pi*alpha*n + theta)
        H[n, n] = V
        if n > 0: H[n, n-1] = 1.0
        if n < N-1: H[n, n+1] = 1.0
    return np.sort(np.linalg.eigvalsh(H))

def log_periodic_analysis(eigs, nbins=128):
    """FFT of log-gap distribution. Returns (freqs, amps)."""
    gaps = np.diff(eigs)
    gaps = gaps[gaps > 1e-8]
    if len(gaps) < 20:
        return None, None
    log_gaps = np.log(gaps)
    hist, edges = np.histogram(log_gaps, bins=nbins, density=True)
    centers = 0.5*(edges[1:] + edges[:-1])
    dlog = centers[1] - centers[0]
    # Remove mean
    dc = hist - hist.mean()
    fft = np.abs(np.fft.rfft(dc))
    freqs = np.fft.rfftfreq(len(hist), d=dlog)
    return freqs, fft

print(f"{'lambda':>8} {'top-3 peak freqs':>30} {'amp':>15}")
print("-"*72)

for lam in [0.5, 1.0, 1.5, 2.0]:
    eigs = aubry_spectrum(N=800, alpha=INV_PHI, lam=lam)
    freqs, fft = log_periodic_analysis(eigs, nbins=128)
    if freqs is None:
        print(f"{lam:>8.2f} no gaps")
        continue
    # Top-3 peaks (excluding DC)
    idx = np.argsort(fft[1:])[::-1][:3] + 1
    peaks = [(freqs[i], fft[i]) for i in idx]
    # Print
    s = "  ".join(f"{f:.3f}" for f, _ in peaks)
    a = peaks[0][1]
    print(f"{lam:>8.2f} {s:>30} {a:>15.4f}")

print()
print("="*72)
print("Сравнение с f_target = 1/log(phi)")
print("="*72)
print(f"f_target = {1/LAM_PHI:.4f}")
print()
print("Если пики FFT близки к f_target или его кратным:")
print("  -> есть phi-log-periodic структура")
print()

# Check specifically at critical point lam=1
eigs = aubry_spectrum(N=800, alpha=INV_PHI, lam=1.0)
freqs, fft = log_periodic_analysis(eigs, nbins=256)
if freqs is not None:
    idx = np.argsort(fft[1:])[::-1][:5] + 1
    print("Top-5 peaks at lam=1.0 (critical):")
    for i in idx:
        ratio = freqs[i] / (1/LAM_PHI) if freqs[i] > 0 else 0
        print(f"  freq = {freqs[i]:.4f}, amp = {fft[i]:.4f}, freq/f_target = {ratio:.4f}")

# Also try alpha=1/phi vs alpha=sqrt(2)-1 (both quadratic irrationals)
print()
print("="*72)
print("Comparison: alpha=1/phi vs alpha=sqrt(2)-1")
print("="*72)
for alpha_name, alpha in [("1/phi", 1/PHI), ("sqrt(2)-1", math.sqrt(2)-1)]:
    eigs = aubry_spectrum(N=800, alpha=alpha, lam=1.0)
    freqs, fft = log_periodic_analysis(eigs, nbins=128)
    if freqs is not None:
        idx = np.argsort(fft[1:])[::-1][:3] + 1
        peaks = [freqs[i] for i in idx]
        print(f"alpha={alpha_name:>12s} ({alpha:.6f}): top peaks = {[f'{p:.3f}' for p in peaks]}")
