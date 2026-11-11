"""Drude–Lorentz optical response models in energy (eV) units.

The model is linear, local, isotropic and causal for non-negative oscillator
strengths and positive damping. It is not a first-principles electronic-structure
calculation and does not include anisotropy or spatial dispersion.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from scipy.integrate import trapezoid

PLANCK_J_S = 6.62607015e-34
HBAR_J_S = PLANCK_J_S / (2.0 * np.pi)
ELEMENTARY_CHARGE_C = 1.602176634e-19
VACUUM_PERMITTIVITY_F_M = 8.8541878128e-12
HC_EV_M = PLANCK_J_S * 299_792_458.0 / ELEMENTARY_CHARGE_C


@dataclass(frozen=True)
class LorentzOscillator:
    """One bound-electron resonance in the Drude–Lorentz model."""

    resonance_energy_ev: float
    damping_ev: float
    strength: float = 1.0


def _validate_energy_grid(energy_ev: np.ndarray) -> np.ndarray:
    energy = np.asarray(energy_ev, dtype=float)
    if energy.ndim != 1 or energy.size < 10:
        raise ValueError("Energy grid must be a one-dimensional array with at least 10 points.")
    if not np.all(np.isfinite(energy)) or np.any(energy <= 0.0):
        raise ValueError("Photon energies must be finite and strictly positive.")
    if np.any(np.diff(energy) <= 0.0):
        raise ValueError("Photon energies must be strictly increasing.")
    return energy


def drude_lorentz_dielectric(
    energy_ev: np.ndarray,
    eps_inf: float,
    plasma_energy_ev: float,
    drude_damping_ev: float,
    oscillators: Sequence[LorentzOscillator] = (),
) -> np.ndarray:
    """Return complex relative dielectric function epsilon1 + i epsilon2.

    In photon-energy units E = hbar*omega:

        epsilon(E) = epsilon_inf
          - Ep^2 / (E^2 + i*Gamma_D*E)
          + sum_j f_j*Ep^2 / (E0_j^2 - E^2 - i*Gamma_j*E)

    A passive medium under the exp(-i omega t) convention has epsilon2 >= 0
    for positive photon energy when all oscillator strengths and damping
    constants are non-negative/positive.
    """
    energy = _validate_energy_grid(energy_ev)
    if not np.isfinite(eps_inf) or eps_inf <= 0.0:
        raise ValueError("epsilon_inf must be finite and positive.")
    if not np.isfinite(plasma_energy_ev) or plasma_energy_ev < 0.0:
        raise ValueError("Plasma energy must be finite and non-negative.")
    if not np.isfinite(drude_damping_ev) or drude_damping_ev <= 0.0:
        raise ValueError("Drude damping must be finite and positive.")

    epsilon = np.full(energy.shape, complex(eps_inf, 0.0), dtype=complex)

    if plasma_energy_ev > 0.0:
        epsilon -= plasma_energy_ev**2 / (
            energy**2 + 1j * drude_damping_ev * energy
        )

    for oscillator in oscillators:
        e0 = float(oscillator.resonance_energy_ev)
        gamma = float(oscillator.damping_ev)
        strength = float(oscillator.strength)
        if not np.all(np.isfinite([e0, gamma, strength])):
            raise ValueError("Oscillator parameters must be finite.")
        if e0 <= 0.0 or gamma <= 0.0 or strength < 0.0:
            raise ValueError(
                "Resonance energy and damping must be positive; strength must be non-negative."
            )
        epsilon += (
            strength * plasma_energy_ev**2
            / (e0**2 - energy**2 - 1j * gamma * energy)
        )

    if not np.all(np.isfinite(epsilon)):
        raise FloatingPointError("Non-finite dielectric response produced; inspect model parameters.")
    return epsilon


def optical_constants(
    energy_ev: np.ndarray,
    epsilon: np.ndarray,
) -> dict[str, np.ndarray]:
    """Calculate complex refractive index, wavelength and absorption coefficient."""
    energy = _validate_energy_grid(energy_ev)
    dielectric = np.asarray(epsilon, dtype=complex)
    if dielectric.shape != energy.shape or not np.all(np.isfinite(dielectric)):
        raise ValueError("Dielectric array must match the energy grid and be finite.")

    refractive_index = np.sqrt(dielectric)
    # Select the passive-medium branch: n >= 0 and k >= 0.
    flip = np.imag(refractive_index) < 0.0
    refractive_index[flip] *= -1.0
    n_real = np.real(refractive_index)
    k_extinction = np.imag(refractive_index)
    wavelength_m = HC_EV_M / energy
    alpha_per_m = 4.0 * np.pi * k_extinction / wavelength_m

    return {
        "n": n_real,
        "k": k_extinction,
        "wavelength_m": wavelength_m,
        "absorption_per_m": alpha_per_m,
    }


def normal_incidence_slab_rt(
    energy_ev: np.ndarray,
    complex_index: np.ndarray,
    thickness_m: float,
    incident_index: float = 1.0,
    substrate_index: float = 1.0,
) -> dict[str, np.ndarray]:
    """Coherent normal-incidence R/T/A of one absorbing film on a substrate.

    The incident and substrate media are assumed lossless and non-magnetic.
    The film is plane-parallel, homogeneous and isotropic. The substrate is
    treated as semi-infinite; its back surface is not included.
    """
    energy = _validate_energy_grid(energy_ev)
    index = np.asarray(complex_index, dtype=complex)
    if index.shape != energy.shape or not np.all(np.isfinite(index)):
        raise ValueError("Complex refractive index must match the energy grid.")
    if not np.isfinite(thickness_m) or thickness_m < 0.0:
        raise ValueError("Film thickness must be finite and non-negative.")
    if incident_index <= 0.0 or substrate_index <= 0.0:
        raise ValueError("Incident and substrate refractive indices must be positive.")
    if np.any(np.imag(index) < -1e-12):
        raise ValueError("Passive film refractive index must have non-negative extinction.")

    wavelength = HC_EV_M / energy
    phase = np.exp(1j * 2.0 * np.pi * index * thickness_m / wavelength)
    phase_round_trip = phase**2

    r01 = (incident_index - index) / (incident_index + index)
    r12 = (index - substrate_index) / (index + substrate_index)
    t01 = 2.0 * incident_index / (incident_index + index)
    t12 = 2.0 * index / (index + substrate_index)

    denominator = 1.0 + r01 * r12 * phase_round_trip
    r_amplitude = (r01 + r12 * phase_round_trip) / denominator
    t_amplitude = t01 * t12 * phase / denominator

    reflectance = np.abs(r_amplitude) ** 2
    transmittance = (substrate_index / incident_index) * np.abs(t_amplitude) ** 2
    absorptance = 1.0 - reflectance - transmittance

    if np.min(absorptance) < -1e-7:
        raise FloatingPointError(
            "Slab energy balance became negative; inspect passive optical parameters."
        )
    absorptance = np.clip(absorptance, 0.0, 1.0)
    return {
        "reflectance": np.real(reflectance),
        "transmittance": np.real(transmittance),
        "absorptance": np.real(absorptance),
    }


def optical_conductivity_s_m(
    energy_ev: np.ndarray,
    epsilon: np.ndarray,
    eps_inf: float,
) -> np.ndarray:
    """Complex optical conductivity sigma = -i*epsilon0*omega*(epsilon-eps_inf)."""
    energy = _validate_energy_grid(energy_ev)
    dielectric = np.asarray(epsilon, dtype=complex)
    if dielectric.shape != energy.shape or not np.all(np.isfinite(dielectric)):
        raise ValueError("Dielectric array must match the energy grid and be finite.")
    if not np.isfinite(eps_inf):
        raise ValueError("epsilon_inf must be finite.")
    omega_rad_s = energy * ELEMENTARY_CHARGE_C / HBAR_J_S
    return (
        -1j
        * VACUUM_PERMITTIVITY_F_M
        * omega_rad_s
        * (dielectric - eps_inf)
    )


def kramers_kronig_eps1_finite_window(
    energy_ev: np.ndarray,
    epsilon2: np.ndarray,
    eps_inf: float,
) -> np.ndarray:
    """Numerically estimate epsilon1 from epsilon2 using a finite-window KK integral.

    Implements the positive-frequency relation

        eps1(E) - eps_inf = (1/pi) PV integral
          [eps2(sqrt(x)) / (x - E^2)] dx,  x = E'^2.

    The singularity is removed by subtracting the local epsilon2 value before
    integration, then adding the analytic finite-boundary logarithm. The first
    and last grid points are returned as NaN because the principal-value
    logarithm is undefined at the truncation boundaries. Spectral tails outside
    the user-selected energy interval are omitted, so this is a diagnostic, not
    an exact reconstruction over an infinite frequency domain.
    """
    energy = _validate_energy_grid(energy_ev)
    imag = np.asarray(epsilon2, dtype=float)
    if imag.shape != energy.shape or not np.all(np.isfinite(imag)):
        raise ValueError("epsilon2 must be finite and match the energy grid.")
    if not np.isfinite(eps_inf):
        raise ValueError("epsilon_inf must be finite.")

    x = energy**2
    derivative = np.gradient(imag, x)
    reconstructed = np.full_like(energy, np.nan, dtype=float)

    for i in range(1, energy.size - 1):
        delta = x - x[i]
        integrand = np.empty_like(imag)
        non_singular = delta != 0.0
        integrand[non_singular] = (imag[non_singular] - imag[i]) / delta[non_singular]
        integrand[i] = derivative[i]

        boundary_term = imag[i] * np.log(
            (x[-1] - x[i]) / (x[i] - x[0])
        )
        reconstructed[i] = eps_inf + (
            trapezoid(integrand, x) + boundary_term
        ) / np.pi

    return reconstructed
