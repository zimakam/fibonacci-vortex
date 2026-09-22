"""Generalized Burgers flow for FibonacciVortex.

∂_t Γ = ν(Γ''/r − Γ'/r²) + a(r)Γ'

с  a/ν = (S₂/S₁)·λ/r + 1/r²,  Sₙ = Σ kⁿ F_k e^{−kλr},  k = 0..4.

BC: Dirichlet с обеих сторон (обе границы удерживаются на Γ_fib).
Оператор согласован с dynamics/eigen_fix.py (интерьер 78 точек,
N=80, r ∈ [0.5, 6.0]).
"""
import numpy as np
import math

PHI = (1.0 + math.sqrt(5.0)) / 2.0
LAM = math.log(PHI)
FIB = [1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144]
K = 5

__all__ = ["BurgersFlow", "a_over_nu", "target_profile", "interior_gap"]


def a_over_nu(r: float) -> float:
    num = sum(k * k * FIB[k] * math.exp(-k * LAM * r) for k in range(K))
    den = sum(k * FIB[k] * math.exp(-k * LAM * r) for k in range(K))
    if abs(den) < 1e-15:
        return 0.0
    return (num / den) * LAM / r + 1.0 / (r * r)


def target_profile(r: np.ndarray) -> np.ndarray:
    return np.array([
        sum(FIB[k] * math.exp(-k * LAM * rr) for k in range(K))
        for rr in r
    ])


class BurgersFlow:
    """Backward-Euler эволюция обобщённого уравнения Бюргера."""

    def __init__(self, N: int = 80, r_min: float = 0.5,
                 r_max: float = 6.0, nu: float = 1.0):
        self.N = int(N)
        self.r = np.linspace(r_min, r_max, self.N)
        self.dr = float(self.r[1] - self.r[0])
        self.nu = float(nu)
        self.G_target = target_profile(self.r)
        self.G = self.G_target.copy()
        self._build_operator()
        self.history = []

    def _build_operator(self) -> None:
        """L такой, что ∂_t Γ = ν · L · Γ."""
        N, dr, r = self.N, self.dr, self.r
        ani = np.array([a_over_nu(rr) for rr in r])
        L = np.zeros((N, N))
        # Граничные строки: Γ_0 и Γ_{N-1} удерживаются (Dirichlet).
        L[0, 0] = 1.0
        L[N - 1, N - 1] = 1.0
        # Интерьер 1..N-2 — та же формула, что в eigen_fix.py.
        for i in range(1, N - 1):
            L[i, i - 1] = (1.0 / (dr * dr * r[i])
                           + 1.0 / (2 * dr * r[i] * r[i])
                           - ani[i] / (2 * dr))
            L[i, i]     = -2.0 / (dr * dr * r[i])
            L[i, i + 1] = (1.0 / (dr * dr * r[i])
                           - 1.0 / (2 * dr * r[i] * r[i])
                           + ani[i] / (2 * dr))
        self.L = L
        self.ani = ani

    def evolve(self, T: float = 0.5, dt: float = 0.005) -> np.ndarray:
        """Неявная эволюция. Границы удерживаются на target."""
        I = np.eye(self.N)
        A = I - dt * self.nu * self.L
        steps = max(1, int(round(T / dt)))
        for _ in range(steps):
            G_new = np.linalg.solve(A, self.G)
            # Dirichlet: принудительно target на границах.
            G_new[0] = self.G_target[0]
            G_new[-1] = self.G_target[-1]
            self.G = G_new
            self.history.append(self.energy())
        return self.G

    def energy(self) -> float:
        delta = self.G - self.G_target
        return float(np.trapezoid(delta * delta * self.r, self.r))

    def interior_gap(self) -> float:
        """−max Re(eig(L_interior)). Должно совпасть с eigen_fix.py."""
        M = self.N - 2
        L_int = self.L[1:self.N - 1, 1:self.N - 1]
        eigs = np.linalg.eigvals(L_int)
        return -float(np.max(np.real(eigs)))

    def summary(self) -> str:
        return (f"BurgersFlow(N={self.N}, r∈[{self.r[0]:.2f},"
                f"{self.r[-1]:.2f}], nu={self.nu}, "
                f"gap={self.interior_gap():.6e})")


def interior_gap(N: int = 80) -> float:
    """Утилита: gap без создания долгой эволюции."""
    bf = BurgersFlow(N=N)
    return bf.interior_gap()


if __name__ == "__main__":
    print("=" * 60)
    print("BURGERS FLOW — VERIFICATION")
    print("=" * 60)
    bf = BurgersFlow(N=80)
    gap = bf.interior_gap()
    print(f"interior gap   : {gap:.6e}")
    print(f"target vs eigen_fix expected: 1.306777e-01")
    print(f"match          : {abs(gap - 1.306777e-01) < 1e-5}")
    E0 = bf.energy()
    bf.evolve(T=0.5, dt=0.005)
    E1 = bf.energy()
    # target стационарен; после эволюции остаётся дискретный
    # остаток O(dr^2), экспоненциально малый — не «рост энергии».
    print(f"E(0)             : {E0:.6e}   (target stationary)")
    print(f"E(0.5)           : {E1:.6e}   (discrete residual)")
    print(f"residual small   : {E1 < 1e-3}")
    print(f"residual bounded : {E1 < 1.0}   (no blow-up)")
