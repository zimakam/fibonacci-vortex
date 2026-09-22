import math
import cmath
import numpy as np

PHI = (1+math.sqrt(5))/2
LAM = math.log(PHI)

# ============================================================
# Riemann-Siegel Z-function
# Z(t) = 2 * sum_{n=1}^{N} cos(theta(t) - t*log(n)) / sqrt(n)
# N = floor(sqrt(t/(2*pi)))
# theta(t) = Im(log Gamma(1/4 + i t/2)) - (t/2)*log(pi)
# ============================================================

def theta(t):
    """Riemann-Siegel theta function."""
    # Use asymptotic formula for large t
    if t < 1:
        return 0.0
    # theta(t) = (t/2)*log(t/(2*pi)) - t/2 - pi/8 + 1/(48t) + ...
    return (t/2)*math.log(t/(2*math.pi)) - t/2 - math.pi/8 + 1/(48*t)

def Z(t):
    """Riemann-Siegel Z-function (real-valued)."""
    if t < 1:
        return 0.0
    N = int(math.sqrt(t/(2*math.pi)))
    if N < 1:
        return 0.0
    th = theta(t)
    s = 0.0
    for n in range(1, N+1):
        s += math.cos(th - t*math.log(n)) / math.sqrt(n)
    return 2.0 * s

def find_zeros(t_min=10.0, t_max=200.0, dt=0.05):
    """Find zeros of Z(t) by sign change."""
    zeros = []
    t = t_min
    prev = Z(t)
    while t < t_max:
        t_next = t + dt
        cur = Z(t_next)
        if prev * cur < 0:
            # Bisection
            a, b = t, t_next
            fa = prev
            for _ in range(60):
                m = 0.5*(a+b)
                fm = Z(m)
                if abs(fm) < 1e-10:
                    break
                if fa * fm < 0:
                    b = m
                else:
                    a = m
                    fa = fm
            zeros.append(0.5*(a+b))
        t = t_next
        prev = cur
    return zeros

print("="*72)
print("RIEMANN: phi-structure in zeta zeros (empirical test)")
print("="*72)
print()

# Known first zeros (Odlyzko reference)
known = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
         37.586178, 40.918719, 43.327073, 48.005151, 49.773832]

print("Computing Z(t) zeros...")
zeros = find_zeros(t_min=10.0, t_max=150.0, dt=0.05)
print(f"Found {len(zeros)} zeros in [10, 150]")
print()

# Compare first 5 with known
print("First 5 zeros (comparison with Odlyzko):")
print(f"{'index':>5} {'found':>12} {'known':>12} {'error':>12}")
print("-"*48)
for i in range(min(5, len(zeros))):
    err = abs(zeros[i] - known[i]) if i < len(known) else 0
    k = known[i] if i < len(known) else 0
    print(f"{i+1:>5} {zeros[i]:>12.6f} {k:>12.6f} {err:>12.2e}")

print()

# Gaps between consecutive zeros
if len(zeros) > 2:
    gaps = np.diff(zeros)
    # Normalize by mean gap
    mean_gap = gaps.mean()
    gaps_norm = gaps / mean_gap

    print("Gap statistics:")
    print(f"  Number of gaps: {len(gaps)}")
    print(f"  Mean gap: {mean_gap:.6f}")
    print(f"  Std gap:  {gaps.std():.6f}")
    print(f"  Min gap:  {gaps.min():.6f}")
    print(f"  Max gap:  {gaps.max():.6f}")
    print()

    # GUE prediction: gaps ~ Wigner surmise
    # P(s) = (32/pi^2) s^2 exp(-4 s^2/pi)
    # Test: fraction of small gaps (< 0.5)
    small = np.sum(gaps_norm < 0.5) / len(gaps_norm)
    print(f"Fraction with normalized gap < 0.5: {small:.4f}")
    print(f"  (GUE predicts ~0.12, Poisson ~0.39)")
    print()

    # Log-periodic analysis: FFT of (gaps - mean) over index
    # If φ-structure: peaks at frequency related to log(φ)
    dg = gaps - gaps.mean()
    if len(dg) > 10:
        fft = np.abs(np.fft.rfft(dg))
        freqs = np.fft.rfftfreq(len(dg), d=1.0)  # in 1/index
        top = np.argsort(fft)[::-1][:5]
        print("Top-5 FFT peaks of (gaps - mean):")
        for i in top:
            print(f"  freq = {freqs[i]:.4f}, amp = {fft[i]:.4f}")
        print()

    # Direct test: are gap ratios close to phi?
    ratios = gaps_norm[1:] / gaps_norm[:-1]
    phi_close = np.sum(np.abs(ratios - PHI) < 0.05)
    print(f"Gap ratios close to phi (within 5%): {phi_close} / {len(ratios)}")
    inv_phi_close = np.sum(np.abs(ratios - 1/PHI) < 0.05)
    print(f"Gap ratios close to 1/phi (within 5%): {inv_phi_close} / {len(ratios)}")
    print()

print("="*72)
print("CONCLUSION")
print("="*72)
print()
print("If GUE stats match and no phi peaks -> zeros are GUE, not phi")
print("If phi peaks found -> potential phi-structure (very interesting)")
