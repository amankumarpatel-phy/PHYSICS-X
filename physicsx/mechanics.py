"""Small, unit-aware mechanics models used by the PHYSICS-X interface.

All inputs are SI values unless the calling interface explicitly converts them.
"""

from __future__ import annotations

import numpy as np


def axial_response(
    force_n: float,
    area_m2: float,
    length_m: float,
    youngs_modulus_pa: float,
    yield_strength_pa: float,
) -> dict[str, float]:
    """Return elementary axial bar response under uniform uniaxial loading.

    Assumptions: prismatic bar, uniform cross-section, small deformation,
    linear-elastic material response, and quasi-static axial load.
    """
    values = (force_n, area_m2, length_m, youngs_modulus_pa, yield_strength_pa)
    if not all(np.isfinite(v) for v in values):
        raise ValueError("All inputs must be finite numbers.")
    if area_m2 <= 0 or length_m <= 0 or youngs_modulus_pa <= 0:
        raise ValueError("Area, length, and Young's modulus must be positive.")
    if yield_strength_pa <= 0:
        raise ValueError("Yield strength must be positive.")

    stress_pa = force_n / area_m2
    strain = stress_pa / youngs_modulus_pa
    elongation_m = strain * length_m
    fos = yield_strength_pa / abs(stress_pa) if stress_pa != 0 else float("inf")

    return {
        "stress_pa": float(stress_pa),
        "strain": float(strain),
        "elongation_m": float(elongation_m),
        "factor_of_safety": float(fos),
        "yield_utilization": float(abs(stress_pa) / yield_strength_pa),
    }


def idealized_bilinear_curve(
    youngs_modulus_pa: float,
    yield_strength_pa: float,
    ultimate_strength_pa: float,
    ultimate_strain: float,
    points: int = 300,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Generate an illustrative engineering stress–strain curve.

    The model is linear elastic up to yield, then linearly strain-hardens to
    the user-specified ultimate point. It does not model necking, damage,
    unloading, rate effects, or a material-specific constitutive law.
    """
    if youngs_modulus_pa <= 0 or yield_strength_pa <= 0:
        raise ValueError("Young's modulus and yield strength must be positive.")
    if ultimate_strength_pa < yield_strength_pa:
        raise ValueError("Ultimate strength must be at least the yield strength.")
    if ultimate_strain <= 0 or points < 10:
        raise ValueError("Ultimate strain must be positive and points must be >= 10.")

    yield_strain = yield_strength_pa / youngs_modulus_pa
    if ultimate_strain <= yield_strain:
        raise ValueError("Ultimate strain must exceed the elastic yield strain.")

    strain = np.linspace(0.0, ultimate_strain, points)
    stress = np.where(
        strain <= yield_strain,
        youngs_modulus_pa * strain,
        yield_strength_pa
        + (ultimate_strength_pa - yield_strength_pa)
        * (strain - yield_strain)
        / (ultimate_strain - yield_strain),
    )
    return strain, stress, float(yield_strain)
