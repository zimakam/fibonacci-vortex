# fibonacci-vortex — research program

Author: Зиявутдинов Магомед Камалович (Zimaka)
ORCID: 0009-0005-9212-9921 · Email: zimakam@gmail.com · License: MIT

Core object: a vortex whose circulation is a Fibonacci-weighted, phi-modulated
hierarchy of K layers:

    Gamma(r) = Gamma_0 * sum_{k=0}^{K-1} F_k * phi^(-k * r / r_c)

Self-similar: Gamma(phi*r) / Gamma(r) -> 1/phi. At K=5 the core is
exactly phi-self-similar (from the c_inf expansion, see proofs/).

This repo is a **research program** — 8 independent lines, ~31k LOC.
It is NOT part of the blockchain project (`zimakam/master-`).
They share only the Gamma(r) formula and the constant c_inf.

---

## Layout

| Folder | Lines | Topic |
|---|---|---|
| millennium/yang_mills/ | ~20 300 | lattice SU(2): quaternion plaquettes, metropolis, wilson_loop |
| millennium/ns/ | ~14 400 | Navier-Stokes: implicit, nonstat, energy, BKM, spectral |
| millennium/pvsnp/ | ~7 100 | 3-SAT DPLL + Fibonacci branching heuristic |
| millennium/riemann/ | ~4 200 | Riemann-Siegel Z(t), theta(t) |
| sparc_fit/ | ~19 000 | fitting 171 SPARC galaxies (vortex vs NFW) |
| dynamics/ | ~13 000 | Burgers flow, eigensolvers (interior 78 pts, N=80) |
| phi_systems/ | ~10 000 | Aubry-Andre localisation (2000x2000 dense) |
| proofs/ | ~4 500 | c_inf = -log(phi) / (2*(phi+2)) — analytic + numeric |
| paper/ | — | LaTeX draft (JPhysA-126040 submitted, declined at desk) |
| fibonacci_vortex.py, theory.py | — | minimal reference implementation |

---

## The constant c_inf

`proofs/c_inf_verify.py` derives, three independent ways:

    c_inf = -log(phi) / (2*(phi + 2)) = -0.066501838644...

Verified to machine precision against:
    A_inf = 1/(3 - phi) = 1/(1 + phi^-2)
    B_inf = -phi^-2 / (1 + phi^-2)^2
via phi^2 = phi + 1  =>  phi^2 + 1 = phi + 2.

The selection K* = 5.43 -> K = 5 follows from alpha(K) = 1
(exact self-similarity condition).

---

## Status of each line

| Line | Code | Status |
|---|---|---|
| FibonacciVortex + c_inf | fibonacci_vortex.py, proofs/ | verified, closed form |
| Yang-Mills | ym_phi.py | runs on small lattices (SU(2), metropolis) |
| Navier-Stokes | ns/step1..4 | BKM integral finite, spectral gap gamma > 0 |
| PvsNP (SAT) | pvsnp_phi.py | Fibonacci branching heuristic NEUTRAL (ratio 0.99) |
| Riemann | riemann_phi.py | Z(t) + theta(t), zeros not proven |
| SPARC | sparc_fit/run_fit*.py | VORTEX vs NFW: wins 134/171; optimiser caveat noted |
| Burgers | dynamics/burgers_flow.py | BC = Gamma_fib, interior gap matches to 7 digits |
| Aubry-Andre | phi_systems/aubry_v2.py | spectrum, log-DOS — no claim yet |

Negative results are recorded honestly (see git log for PvsNP).

---

## Usage

    python3 fibonacci_vortex.py --selftest
    python3 theory.py
    python3 proofs/c_inf_verify.py
    python3 millennium/yang_mills/ym_phi.py
    python3 dynamics/burgers_flow.py

Dependencies: numpy, math (stdlib). Some scripts use matplotlib (Agg backend).

---

## Related

- `zimakam/fibonacci-vortex-core` — minimal, published core (MIT)
- `zimakam/master-` — blockchain prototype (uses only the Gamma(r) formula)
- `zimakam/pqpchain-` (private) — PQC toolkit

## Not part of the blockchain

The blockchain project uses exactly two objects from here:
    Gamma(r)       — as a source of channel speeds (vortex_consensus.py)
    c_inf          — as a regression constant (theory_vortex.py)
Everything else (millennium, SPARC, Burgers, Aubry) is independent.
