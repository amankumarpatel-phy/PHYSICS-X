"""Harmonic phonons of a one-dimensional diatomic chain.

The model contains two masses per primitive cell and equal nearest-neighbour
harmonic springs. It is a pedagogical 1D model, not a 3D material predictor.
"""

from __future__ import annotations

import numpy as np

ATOMIC_MASS_UNIT_KG = 1.66053906660e-27


def diatomic_chain_dispersion(
    lattice_constant_m: float,
    force_constant_n_m: float,
    mass1_kg: float,
    mass2_kg: float,
    q_points: int = 801,
) -> dict[str, np.ndarray]:
    """Calculate acoustic and optical branches across the 1D first BZ.

    Parameters use SI units. lattice_constant_m is the primitive-cell repeat
    length (the two sites are separated by a/2 in this model).
    The returned eigenfrequencies are angular frequencies in rad/s.
    """
    values = (lattice_constant_m, force_constant_n_m, mass1_kg, mass2_kg)
    if not all(np.isfinite(value) for value in values):
        raise ValueError("Physical parameters must be finite.")
    if any(value <= 0.0 for value in values):
        raise ValueError("Lattice constant, force constant and masses must be positive.")
    if not isinstance(q_points, int) or q_points < 3:
        raise ValueError("q_points must be an integer of at least 3.")

    q = np.linspace(-np.pi / lattice_constant_m, np.pi / lattice_constant_m, q_points)
    acoustic = np.empty_like(q)
    optical = np.empty_like(q)
    coupling_scale = np.sqrt(mass1_kg * mass2_kg)

    for index, wavevector in enumerate(q):
        phase = np.exp(-1j * wavevector * lattice_constant_m)
        off_diagonal = -force_constant_n_m * (1.0 + phase) / coupling_scale
        dynamical_matrix = np.array(
            [
                [2.0 * force_constant_n_m / mass1_kg, off_diagonal],
                [np.conjugate(off_diagonal), 2.0 * force_constant_n_m / mass2_kg],
            ],
            dtype=complex,
        )
        eigenvalues = np.linalg.eigvalsh(dynamical_matrix)
        # The translational acoustic mode is exactly zero at Gamma. Enforce
        # that analytical result to avoid loss of precision in the determinant.
        if wavevector == 0.0:
            eigenvalues[0] = 0.0
        eigenvalues = np.maximum(eigenvalues, 0.0)
        frequencies = np.sqrt(eigenvalues)
        acoustic[index], optical[index] = frequencies

    return {
        "q_rad_m": q,
        "q_reduced": q * lattice_constant_m / np.pi,
        "acoustic_omega_rad_s": acoustic,
        "optical_omega_rad_s": optical,
    }


def angular_frequency_to_thz(omega_rad_s: np.ndarray) -> np.ndarray:
    """Convert angular frequency in rad/s to ordinary frequency in THz."""
    omega = np.asarray(omega_rad_s, dtype=float)
    if not np.all(np.isfinite(omega)) or np.any(omega < 0.0):
        raise ValueError("Angular frequencies must be finite and non-negative.")
    return omega / (2.0 * np.pi * 1.0e12)


def phonon_density_of_states(
    acoustic_omega_rad_s: np.ndarray,
    optical_omega_rad_s: np.ndarray,
    bins: int = 100,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Estimate 1D phonon DOS by histogramming sampled q-point frequencies.

    The returned DOS is in states / THz / primitive cell. Its integral over
    frequency is 2, corresponding to two vibrational branches per cell.
    """
    acoustic = np.asarray(acoustic_omega_rad_s, dtype=float)
    optical = np.asarray(optical_omega_rad_s, dtype=float)
    if acoustic.ndim != 1 or optical.ndim != 1 or acoustic.shape != optical.shape:
        raise ValueError("Acoustic and optical branches must be equal-length 1D arrays.")
    if acoustic.size < 2:
        raise ValueError("At least two q-point samples are required.")
    if not np.all(np.isfinite(acoustic)) or not np.all(np.isfinite(optical)):
        raise ValueError("Frequencies must be finite.")
    if np.any(acoustic < 0.0) or np.any(optical < 0.0):
        raise ValueError("Frequencies cannot be negative.")
    if not isinstance(bins, int) or bins < 10:
        raise ValueError("bins must be an integer of at least 10.")

    frequencies = np.concatenate(
        [angular_frequency_to_thz(acoustic), angular_frequency_to_thz(optical)]
    )
    maximum = float(np.max(frequencies))
    if maximum <= 0.0:
        raise ValueError("At least one phonon frequency must be positive.")

    edges = np.linspace(0.0, maximum * 1.000001, bins + 1)
    bin_width = edges[1] - edges[0]
    # Each branch contributes one state per sampled primitive-cell wavevector.
    weights = np.full(frequencies.size, 1.0 / acoustic.size)
    counts, edges = np.histogram(frequencies, bins=edges, weights=weights)
    dos = counts / bin_width
    centres = 0.5 * (edges[:-1] + edges[1:])
    return centres, dos, edges
