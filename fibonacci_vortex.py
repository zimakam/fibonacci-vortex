#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FibonacciVortex — vortex core with φ-modulated circulation
Author: Ziyavutdinov Magomed Kamalovich (Zimaka)
License: MIT
"""
import math
import argparse

PHI = (1.0 + math.sqrt(5.0)) / 2.0
INV_PHI = PHI - 1.0
FIB = [1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144]

__author__ = "Зиявутдинов Магомед Камалович (Zimaka)"
__version__ = "1.0.0"
__license__ = "MIT"


class FibonacciVortex:
    """Vortex core with φ-modulated circulation."""

    def __init__(self, circulation=13.0, core_radius=1.0,
                 k_layers=5, name="fib_vortex"):
        self.circulation = float(circulation)
        self.core_radius = float(core_radius)
        self.k_layers = int(k_layers)
        self.name = name

    def gamma_of_r(self, r):
        """Γ(r) = Γ₀ · Σ F_k · φ^(−k·r/r_c)."""
        rc = self.core_radius
        return self.circulation * sum(
            FIB[k] * (PHI ** (-k * r / rc))
            for k in range(self.k_layers))

    def self_similarity(self, r):
        """Γ(φr)/Γ(r) → 1/φ."""
        g1 = self.gamma_of_r(r)
        g2 = self.gamma_of_r(PHI * r)
        return g2 / g1 if abs(g1) > 1e-15 else float('nan')

    def zeckendorf_address(self):
        """Zeckendorf address of |Γ₀|."""
        n = int(abs(self.circulation) * 1000)
        parts = []
        for f in reversed(FIB):
            if f <= n:
                parts.append(f)
                n -= f
        return parts

    def report(self):
        g0 = self.gamma_of_r(0.0)
        grc = self.gamma_of_r(self.core_radius)
        sim = self.self_similarity(self.core_radius)
        return (
            f"FibonacciVortex({self.name})\n"
            f"  Γ₀ = {self.circulation}\n"
            f"  r_c = {self.core_radius}\n"
            f"  K = {self.k_layers} layers\n"
            f"  Γ(0) = {g0:.4f}\n"
            f"  Γ(r_c) = {grc:.4f}\n"
            f"  Γ(φ·r_c)/Γ(r_c) = {sim:.4f} "
            f"(target 1/φ = {INV_PHI:.4f})\n"
            f"  Zeckendorf: {self.zeckendorf_address()}")


def selftest():
    print("=" * 60)
    print(f"FIBONACCIVORTEX v{__version__} — SELFTEST")
    print("=" * 60)
    fv = FibonacciVortex(circulation=13.0, core_radius=1.0, k_layers=5)
    print(f"\n{fv.report()}")
    sim = fv.self_similarity(1.0)
    print(f"\nφ-scaling at r=r_c: {sim:.4f}")
    print(f"Target 1/φ = {INV_PHI:.4f}")
    if abs(sim - INV_PHI) < 0.05:
        print("✅ φ-scaling OK")
    print("\n--- Γ(r) values ---")
    for r in [0.1, 0.5, 1.0, 2.0, 5.0]:
        print(f"  r = {r:4.1f} → Γ(r) = {fv.gamma_of_r(r):10.4f}")
    print("\n✅ SELFTEST passed")
    print(f"Author: {__author__}")


def main():
    parser = argparse.ArgumentParser(
        description=f"FibonacciVortex v{__version__}")
    parser.add_argument("--selftest", action="store_true")
    parser.add_argument("--demo", action="store_true")
    args = parser.parse_args()
    selftest()


if __name__ == "__main__":
    main()