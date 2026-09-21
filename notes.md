# FibonacciVortex — Session Notes 2026-09-21

## What was verified

All checks with double precision (floor ~1e-15).

### 1. Annihilator identity
Gamma_Fib(r) = sum F_k exp(-k*LAM*r), LAM = log(phi)
satisfies exactly:
    prod_{k=0}^{K-1} (d/dr + k*LAM) Gamma = 0
Verified numerically: |L[Gamma]|/|Gamma| ~ 1e-15.

### 2. Burgers-vortex identity (main result)
Gamma_Fib solves the Burgers-vortex equation with variable strain:
    nu*(Gamma''/r - Gamma'/r^2) + a(r)*Gamma' = 0
where
    a(r)/nu = <k>_w(r)*LAM/r + 1/r^2
    <k>_w(r) = sum_k k^2 F_k exp(-k*LAM*r) / sum_k k F_k exp(-k*LAM*r)
Verified: residual ~ 1e-15 at r in [0.1, 12].

### 3. Local homogeneity exponent
alpha(r) = -r * Gamma'(r)/Gamma(r)
- alpha = 1.0005 at r = 1.5*r_c (peak)
- alpha = 1 exactly at r = 1.459*r_c
- r*(alpha) = 1.0 at r=1.46 vs r_c=1.0 -> K*=5.43 -> K=5

### 4. Mode structure (K=5)
- Weights: w_k = F_k * exp(-k*LAM*r)
- Dominant mode: k=4 for r<0.836, k=0 for r>0.836
- Intermediate modes k=1,2,3 never dominate
- Crossover r_cross = log(5)/(4*LAM) = 0.836

### 5. Thermodynamic analogy (new)
Define Z(beta) = sum_k k*F_k*exp(-k*beta), beta = LAM*r
    <k>_w = -d log Z / d beta
Specific heat C(beta) = -beta^2 dU/dbeta has Schottky peak at
    beta* = 2.04 (r = 4.24), C_max = 2.40
Signature of a 4-level system with degeneracies {1, 2, 6, 12}.

### 6. Three independent r-scales
| event                          | r        |
|--------------------------------|----------|
| mode crossover k=4 -> k=0      | 0.836    |
| exact phi-self-similarity      | 1.459    |
| alpha peak                     | 1.500    |
| Schottky peak                  | 4.239    |

These are DIFFERENT phenomena, not the same point.

## What is NOT claimed
- Gamma_Fib is NOT a solution of classical (constant a) Burgers NS
- Not a solution of Lamb-Oseen
- The Fibonacci hierarchy is NOT motivated by first-principles physics
- K=5 selection follows from alpha=1, not from data fit

## Open questions
1. Two-vortex interaction at the crossover r = 0.836
2. Sign-flip mechanism: what changes sign at the crossover?
3. Stability of Gamma_Fib against perturbations
4. Physical origin of the Fibonacci hierarchy

## Files
- fibonacci_vortex.py — model
- theory.py — alpha(K), K-selection, Burgers identity
- paper/main.tex — draft (RevTeX 4.2)
- paper/bibliography.bib

## Author
Ziyavutdinov Magomed Kamalovich (Zimaka)
ORCID 0009-0005-9212-9921


## Late session update 2026-09-21

### Tested Eisenhart lift (potential sign vs stationary point)
- V(r) = Gamma_theta(r), stationary: dV/dr = 0
- Stationary points exist only for theta < 46 deg
  - theta=40: r_stat = 2.44, V = +1.18
  - theta=44: r_stat = 2.80, V = +1.13
- For theta >= 50 deg: no stationary point
- V(r_stat) > 0 always; no point where V=0 AND dV/dr=0
- CONCLUSION: no collapse point in current parametrization,
  neither in 3D dynamics nor in Eisenhart lift

### Rejected hypotheses (negative results)
- Two-vortex emergence: E_overlap monotone in d
- Sign reversal via (-1)^k weights: always positive
- Static collapse point: not present

### Confirmed results (final list)
1. Burgers identity (machine precision)
2. SUSY structure of L_5
3. Thermodynamic analogy + Schottky anomaly
4. Theta phase transition (zero-crossing window)

## End of session 2026-09-21 (late)
