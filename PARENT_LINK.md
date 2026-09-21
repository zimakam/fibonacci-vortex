# PARENT_LINK.md — FibonacciVortex ↔ Quantum Master v12.0

## What this document is
A one-page bridge between the FibonacciVortex theory in this
repository and the parent system Quantum Master v12.0 (single-file).

## Exact connection points

### 1. Class identity
Parent file contains:
    class FibonacciVortex(VortexCore):
        def gamma_of_r(self, r):
            return self.circulation * sum(
                FIB[k] * (PHI ** (-k * r / rc))
                for k in range(self.k_layers))

This is IDENTICAL to fibonacci_vortex.py in this repo.
The theory in this repo is the mathematical foundation
for that one component of the parent.

### 2. Merkaba identity
Parent: Merkaba(mode='asymmetric')
    v_upper = PHI
    v_lower = INV_PHI
    gamma_upper = 1/sqrt(1 - PHI^2)  -> capped (GAMMA_MAX)
    gamma_lower = 1/sqrt(1 - 1/PHI^2) = sqrt(PHI)

THEOREM (verified 2026-09-21 in this repo):
    gamma(1/phi) = sqrt(phi)   EXACT

### 3. Default parameters
Parent: k_layers=5, FIB = [1,1,2,3,5,8,13,...]
Theory in this repo:
    alpha(K) = (K-1)*log(phi)/2 + c_inf
    alpha = 1 -> K* = 5.43 -> K = 5
Selection of K=5 is derivable, not arbitrary.

### 4. GravitationalVortex
Parent: rho(x) = kappa * sum |omega|^2; nabla^2 phi = 4 pi G rho
Repo: Gamma_Fib solves the generalized Burgers equation:
    nu*(Gamma''/r - Gamma'/r^2) + a(r)*Gamma' = 0
    a/nu = <k>_w*lambda/r + 1/r^2
The 'vortex' in the parent is the physical object;
the Burgers identity in the repo is its mathematical description.

## What this means
- The parent is a computational architecture.
- The repo is a mathematical theory for one of its cores.
- Both were derived independently and converge on phi-structure.
- This is a consistency check, not a coincidence.

## What to do next (priority)
1. Do NOT expand the parent tonight.
2. Finish the FibonacciVortex paper (paper/main.tex).
3. Add citation to parent project in the paper.
4. After first publication: separate paper on Merkaba architecture.

## Status
Bridge documented. No new code. Sleep.
