# fibonacci-vortex
Novel vortex core with φ-modulated circulation. Γ(r) = Γ₀·Σ F_k·φ^(−k·r/r_c). Self-similar with ratio 1/φ. No analogues in literature. MIT.
# FibonacciVortex

> Novel vortex core with φ-modulated circulation.
> Γ(r) = Γ₀ · Σ F_k · φ^(−k·r/r_c) — no analogues in literature.

**Author:** Зиявутдинов Магомед Камалович (Zimaka)
**Email:** zimakam@gmail.com
**ORCID:** 0009-0005-9212-9921
**License:** MIT
**Version:** 1.0.0

---

## What is this

Vortex core model where circulation depends on distance via Fibonacci-weighted φ-exponential:

Γ(r) = Γ₀ · Σ F_k · φ^(−k·r/r_c)

Unlike classical models (Rankine 1865, Lamb-Oseen 1911, Gaussian), where Γ = const, here circulation is self-similar with scale φ.

## Uniqueness

### 1. φ-modulated circulation
Each next layer is φ times weaker, weighted by F_k (Fibonacci).

### 2. Self-similarity
Γ(φr) / Γ(r) → 1/φ

### 3. Zeckendorf addressing
Each vortex has unique Zeckendorf address.

## Quick start

pip install numpy
python fibonacci_vortex.py --selftest

## Author

- Зиявутдинов Магомед Камалович (Zimaka)
- ORCID: 0009-0005-9212-9921
- GitHub: @zimakam

## License

MIT