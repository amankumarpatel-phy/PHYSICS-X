"""Ideal cubic crystal structures and basic geometric properties."""

from __future__ import annotations

import numpy as np

STRUCTURES = {
    "Simple Cubic (SC)": {
        "atoms_per_cell": 1,
        "coordination_number": 6,
        "packing_fraction": np.pi / 6.0,
        "basis": [(0.0, 0.0, 0.0)],
        "contact_factor": 1.0,  # a = 2r
    },
    "Body-Centred Cubic (BCC)": {
        "atoms_per_cell": 2,
        "coordination_number": 8,
        "packing_fraction": np.sqrt(3.0) * np.pi / 8.0,
        "basis": [(0.0, 0.0, 0.0), (0.5, 0.5, 0.5)],
        "contact_factor": 2.0 / np.sqrt(3.0),  # a = 4r/sqrt(3)
    },
    "Face-Centred Cubic (FCC)": {
        "atoms_per_cell": 4,
        "coordination_number": 12,
        "packing_fraction": np.sqrt(2.0) * np.pi / 6.0,
        "basis": [
            (0.0, 0.0, 0.0),
            (0.0, 0.5, 0.5),
            (0.5, 0.0, 0.5),
            (0.5, 0.5, 0.0),
        ],
        "contact_factor": 1.0 / np.sqrt(2.0),  # a = 2sqrt(2)r
    },
}


def unit_cell_points(structure: str, lattice_constant_m: float) -> np.ndarray:
    """Return Cartesian basis-point coordinates in metres."""
    if structure not in STRUCTURES:
        raise ValueError(f"Unsupported structure: {structure}")
    if lattice_constant_m <= 0:
        raise ValueError("Lattice constant must be positive.")
    basis = np.asarray(STRUCTURES[structure]["basis"], dtype=float)
    return basis * lattice_constant_m


def theoretical_density_kg_m3(
    structure: str, lattice_constant_m: float, molar_mass_kg_mol: float
) -> float:
    """Compute ideal crystal density from conventional unit-cell contents."""
    if structure not in STRUCTURES:
        raise ValueError(f"Unsupported structure: {structure}")
    if lattice_constant_m <= 0 or molar_mass_kg_mol <= 0:
        raise ValueError("Lattice constant and molar mass must be positive.")

    avogadro = 6.02214076e23
    atoms = STRUCTURES[structure]["atoms_per_cell"]
    mass_cell = atoms * molar_mass_kg_mol / avogadro
    volume_cell = lattice_constant_m**3
    return float(mass_cell / volume_cell)


def atomic_radius_m(structure: str, lattice_constant_m: float) -> float:
    """Infer hard-sphere atomic radius from the ideal contact geometry."""
    if structure not in STRUCTURES or lattice_constant_m <= 0:
        raise ValueError("Choose a supported structure and positive lattice constant.")
    return float(lattice_constant_m * STRUCTURES[structure]["contact_factor"] / 2.0)
