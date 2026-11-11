"""Advanced crystallography: reciprocal space, Miller planes and powder XRD."""

from __future__ import annotations

from itertools import product

import numpy as np
from scipy.spatial import ConvexHull, HalfspaceIntersection

AVOGADRO = 6.02214076e23


def lattice_geometry(
    structure: str,
    lattice_a_m: float,
    lattice_c_m: float | None = None,
) -> dict[str, np.ndarray | float | str]:
    """Return conventional-cell geometry and atomic basis in fractional units.

    Cubic SC/BCC/FCC cells are conventional cubic cells. HCP uses a primitive
    hexagonal parallelepiped with a two-atom basis.
    """
    if not np.isfinite(lattice_a_m) or lattice_a_m <= 0.0:
        raise ValueError("Lattice parameter a must be finite and positive.")

    if structure in {"Simple Cubic (SC)", "Body-Centred Cubic (BCC)", "Face-Centred Cubic (FCC)"}:
        direct = np.array(
            [[lattice_a_m, 0.0, 0.0],
             [0.0, lattice_a_m, 0.0],
             [0.0, 0.0, lattice_a_m]],
            dtype=float,
        )
        if structure == "Simple Cubic (SC)":
            basis = np.array([[0.0, 0.0, 0.0]])
            atoms_per_cell = 1
        elif structure == "Body-Centred Cubic (BCC)":
            basis = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
            atoms_per_cell = 2
        else:
            basis = np.array(
                [[0.0, 0.0, 0.0],
                 [0.0, 0.5, 0.5],
                 [0.5, 0.0, 0.5],
                 [0.5, 0.5, 0.0]]
            )
            atoms_per_cell = 4
        primitive_volume = abs(float(np.linalg.det(direct))) / atoms_per_cell
    elif structure == "Hexagonal Close-Packed (HCP)":
        c_value = lattice_c_m
        if c_value is None or not np.isfinite(c_value) or c_value <= 0.0:
            raise ValueError("HCP requires a finite, positive c lattice parameter.")
        direct = np.array(
            [[lattice_a_m, 0.0, 0.0],
             [-0.5 * lattice_a_m, np.sqrt(3.0) * 0.5 * lattice_a_m, 0.0],
             [0.0, 0.0, c_value]],
            dtype=float,
        )
        basis = np.array([[0.0, 0.0, 0.0], [2.0 / 3.0, 1.0 / 3.0, 0.5]])
        atoms_per_cell = 2
        primitive_volume = abs(float(np.linalg.det(direct)))
    else:
        raise ValueError(f"Unsupported crystal structure: {structure}")

    reciprocal = 2.0 * np.pi * np.linalg.inv(direct).T
    return {
        "structure": structure,
        "direct": direct,
        "reciprocal": reciprocal,
        "basis": basis,
        "atoms_per_cell": atoms_per_cell,
        "cell_volume_m3": abs(float(np.linalg.det(direct))),
        "primitive_volume_m3": primitive_volume,
    }


def fractional_to_cartesian(fractional: np.ndarray, direct: np.ndarray) -> np.ndarray:
    """Convert row-wise fractional coordinates to Cartesian coordinates."""
    fractional = np.asarray(fractional, dtype=float)
    if fractional.shape[-1] != 3:
        raise ValueError("Fractional coordinates must have a final dimension of 3.")
    return fractional @ np.asarray(direct, dtype=float)


def reciprocal_lattice_points(
    reciprocal: np.ndarray,
    structure: str,
    basis: np.ndarray,
    index_limit: int = 2,
) -> tuple[np.ndarray, np.ndarray]:
    """Generate reciprocal-lattice points, excluding systematic absences."""
    if not isinstance(index_limit, int) or index_limit < 1:
        raise ValueError("index_limit must be an integer of at least 1.")
    indices = np.array(
        [hkl for hkl in product(range(-index_limit, index_limit + 1), repeat=3)
         if hkl != (0, 0, 0)],
        dtype=int,
    )
    if structure != "Hexagonal Close-Packed (HCP)":
        factors = structure_factors(indices, basis)
        indices = indices[np.abs(factors) > 1e-8]
    vectors = indices @ np.asarray(reciprocal, dtype=float)
    return indices, vectors


def structure_factors(hkl: np.ndarray, basis: np.ndarray) -> np.ndarray:
    """Monatomic-basis kinematic structure factor with equal atomic form factors.

    Atomic form factors, Debye-Waller factors, anomalous dispersion and
    polarization corrections are intentionally not included.
    """
    indices = np.asarray(hkl, dtype=int)
    positions = np.asarray(basis, dtype=float)
    phases = 2.0j * np.pi * (indices @ positions.T)
    return np.sum(np.exp(phases), axis=1)


def miller_d_spacing_m(hkl: tuple[int, int, int] | np.ndarray, reciprocal: np.ndarray) -> float:
    """Calculate d_hkl = 2 pi / |G_hkl| using a 2-pi reciprocal convention."""
    indices = np.asarray(hkl, dtype=float)
    if indices.shape != (3,) or not np.all(np.isfinite(indices)) or np.all(indices == 0):
        raise ValueError("Miller indices must be three finite values, not all zero.")
    g_vector = indices @ np.asarray(reciprocal, dtype=float)
    magnitude = float(np.linalg.norm(g_vector))
    if magnitude <= 0.0:
        raise ValueError("Reciprocal vector magnitude must be positive.")
    return 2.0 * np.pi / magnitude


def _cube_edges() -> list[tuple[np.ndarray, np.ndarray]]:
    corners = np.array(list(product((0.0, 1.0), repeat=3)), dtype=float)
    edges = []
    for i, first in enumerate(corners):
        for second in corners[i + 1:]:
            if np.count_nonzero(np.abs(first - second) > 0.0) == 1:
                edges.append((first, second))
    return edges


def miller_plane_section(
    hkl: tuple[int, int, int],
    offset: float,
    direct: np.ndarray,
) -> np.ndarray:
    """Return ordered Cartesian vertices where h*u+k*v+l*w=offset cuts a cell."""
    hkl_array = np.asarray(hkl, dtype=float)
    if hkl_array.shape != (3,) or np.all(hkl_array == 0):
        raise ValueError("Provide three Miller indices, not all zero.")
    if not np.isfinite(offset):
        raise ValueError("Plane offset must be finite.")

    intersections: list[np.ndarray] = []
    tolerance = 1e-10
    for first, second in _cube_edges():
        f1 = float(np.dot(hkl_array, first) - offset)
        f2 = float(np.dot(hkl_array, second) - offset)
        if abs(f1) <= tolerance:
            intersections.append(first)
        if abs(f2) <= tolerance:
            intersections.append(second)
        if f1 * f2 < 0.0:
            fraction = f1 / (f1 - f2)
            intersections.append(first + fraction * (second - first))

    if not intersections:
        raise ValueError("This plane does not intersect the unit cell.")

    unique = np.unique(np.round(np.asarray(intersections), decimals=10), axis=0)
    if len(unique) < 3:
        raise ValueError("The plane intersects the cell in fewer than three points.")

    cartesian = fractional_to_cartesian(unique, direct)
    centre = np.mean(cartesian, axis=0)
    reciprocal_normal = hkl_array @ (2.0 * np.pi * np.linalg.inv(direct).T)
    normal = reciprocal_normal / np.linalg.norm(reciprocal_normal)
    reference = np.array([1.0, 0.0, 0.0])
    if abs(float(np.dot(reference, normal))) > 0.9:
        reference = np.array([0.0, 1.0, 0.0])
    basis_u = np.cross(normal, reference)
    basis_u /= np.linalg.norm(basis_u)
    basis_v = np.cross(normal, basis_u)
    relative = cartesian - centre
    angles = np.arctan2(relative @ basis_v, relative @ basis_u)
    return cartesian[np.argsort(angles)]


def first_brillouin_zone(
    reciprocal: np.ndarray,
    structure: str,
    basis: np.ndarray,
    index_limit: int = 3,
) -> dict[str, np.ndarray | float]:
    """Construct the first Brillouin zone as the reciprocal-lattice Wigner-Seitz cell.

    The half-space for each reciprocal vector G is k.G <= |G|^2/2.
    Centering extinctions are removed for conventional cubic BCC/FCC cells;
    HCP uses the primitive hexagonal reciprocal lattice.
    """
    _, vectors = reciprocal_lattice_points(
        reciprocal, structure, basis, index_limit=index_limit
    )
    if len(vectors) < 6:
        raise ValueError("Not enough reciprocal vectors to construct the zone.")
    halfspaces = np.column_stack(
        [vectors, -0.5 * np.einsum("ij,ij->i", vectors, vectors)]
    )
    intersection = HalfspaceIntersection(halfspaces, np.zeros(3))
    vertices = intersection.intersections
    hull = ConvexHull(vertices)
    return {
        "vertices": vertices,
        "simplices": hull.simplices,
        "volume_m3": float(hull.volume),
        "area_m2": float(hull.area),
    }


def powder_xrd_pattern(
    reciprocal: np.ndarray,
    basis: np.ndarray,
    structure: str,
    wavelength_m: float = 1.5406e-10,
    max_index: int = 6,
    two_theta_min_deg: float = 10.0,
    two_theta_max_deg: float = 120.0,
    fwhm_deg: float = 0.16,
    grid_step_deg: float = 0.02,
) -> tuple[np.ndarray, np.ndarray, list[dict[str, float | str]]]:
    """Build a simple kinematic powder pattern with Gaussian peak broadening.

    Peak weights sum |F_hkl|^2 over symmetry-equivalent points found on the
    sampled reciprocal grid. Equal atomic form factors are assumed.
    """
    if wavelength_m <= 0.0 or not np.isfinite(wavelength_m):
        raise ValueError("X-ray wavelength must be finite and positive.")
    if max_index < 1:
        raise ValueError("max_index must be at least 1.")
    if fwhm_deg <= 0.0 or grid_step_deg <= 0.0:
        raise ValueError("Peak width and grid step must be positive.")
    if two_theta_max_deg <= two_theta_min_deg:
        raise ValueError("2theta maximum must exceed its minimum.")

    raw = np.array(
        [hkl for hkl in product(range(-max_index, max_index + 1), repeat=3)
         if hkl != (0, 0, 0)],
        dtype=int,
    )
    factors = structure_factors(raw, basis)
    g_vectors = raw @ np.asarray(reciprocal, dtype=float)
    g_magnitudes = np.linalg.norm(g_vectors, axis=1)
    d_values = 2.0 * np.pi / g_magnitudes
    sin_theta = wavelength_m / (2.0 * d_values)
    valid = (np.abs(factors) > 1e-8) & (sin_theta <= 1.0)
    raw = raw[valid]
    factors = factors[valid]
    d_values = d_values[valid]
    sin_theta = sin_theta[valid]
    two_theta = np.degrees(2.0 * np.arcsin(sin_theta))
    intensities = np.abs(factors) ** 2

    grouped: dict[float, dict[str, float | str]] = {}
    for hkl, d_value, angle, intensity in zip(raw, d_values, two_theta, intensities):
        if angle < two_theta_min_deg or angle > two_theta_max_deg:
            continue
        key = round(float(angle), 5)
        if key not in grouped:
            grouped[key] = {
                "two_theta_deg": float(angle),
                "d_angstrom": float(d_value * 1e10),
                "intensity": 0.0,
                "representative_hkl": f"({hkl[0]} {hkl[1]} {hkl[2]})",
            }
        grouped[key]["intensity"] = float(grouped[key]["intensity"]) + float(intensity)

    peaks = sorted(grouped.values(), key=lambda item: float(item["two_theta_deg"]))
    x = np.arange(two_theta_min_deg, two_theta_max_deg + grid_step_deg / 2.0, grid_step_deg)
    y = np.zeros_like(x)
    sigma = fwhm_deg / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    for peak in peaks:
        centre = float(peak["two_theta_deg"])
        weight = float(peak["intensity"])
        y += weight * np.exp(-0.5 * ((x - centre) / sigma) ** 2)
    if y.size and float(np.max(y)) > 0.0:
        y = 100.0 * y / np.max(y)
    return x, y, peaks
