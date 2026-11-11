"""Analytic non-relativistic 1D infinite square-well model."""

from __future__ import annotations

import numpy as np

PLANCK = 6.62607015e-34  # J s, exact SI value
HBAR = PLANCK / (2.0 * np.pi)
E_CHARGE = 1.602176634e-19  # C, exact SI value


def energy_levels_ev(width_m: float, mass_kg: float, n_max: int = 6) -> np.ndarray:
    """Return E_n in eV for a well from x=0 to x=L."""
    if width_m <= 0 or mass_kg <= 0 or n_max < 1:
        raise ValueError("Width, mass, and n_max must be positive.")
    n = np.arange(1, n_max + 1, dtype=float)
    energy_j = (n**2 * np.pi**2 * HBAR**2) / (2.0 * mass_kg * width_m**2)
    return energy_j / E_CHARGE


def stationary_wavefunction(
    width_m: float, n: int, x_m: np.ndarray
) -> np.ndarray:
    """Normalized spatial eigenfunction psi_n(x) inside an infinite well.

    The function is zero outside the interval [0, L].
    """
    if width_m <= 0 or n < 1:
        raise ValueError("Width and quantum number n must be positive.")
    x = np.asarray(x_m, dtype=float)
    psi = np.zeros_like(x)
    inside = (x >= 0.0) & (x <= width_m)
    psi[inside] = np.sqrt(2.0 / width_m) * np.sin(
        n * np.pi * x[inside] / width_m
    )
    # Enforce the infinite-wall boundary condition exactly; sin(n*pi) can
    # otherwise leave a small floating-point residual at x = L.
    psi[x == 0.0] = 0.0
    psi[x == width_m] = 0.0
    return psi
