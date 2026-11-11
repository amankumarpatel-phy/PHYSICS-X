import unittest

import numpy as np

from physicsx.crystal import STRUCTURES, atomic_radius_m, theoretical_density_kg_m3, unit_cell_points
from physicsx.mechanics import axial_response, idealized_bilinear_curve
from physicsx.quantum import energy_levels_ev, stationary_wavefunction


class MechanicsTests(unittest.TestCase):
    def test_axial_response(self):
        result = axial_response(
            force_n=1000.0,
            area_m2=1.0e-4,
            length_m=1.0,
            youngs_modulus_pa=200.0e9,
            yield_strength_pa=250.0e6,
        )
        self.assertAlmostEqual(result["stress_pa"], 10.0e6)
        self.assertAlmostEqual(result["strain"], 5.0e-5)
        self.assertAlmostEqual(result["elongation_m"], 5.0e-5)
        self.assertAlmostEqual(result["factor_of_safety"], 25.0)

    def test_zero_load_has_infinite_safety_factor(self):
        result = axial_response(0.0, 1.0e-4, 1.0, 200.0e9, 250.0e6)
        self.assertTrue(np.isinf(result["factor_of_safety"]))

    def test_idealized_curve_reaches_requested_end_point(self):
        strain, stress, yield_strain = idealized_bilinear_curve(
            200.0e9, 250.0e6, 400.0e6, 0.15
        )
        self.assertAlmostEqual(stress[0], 0.0)
        self.assertAlmostEqual(stress[-1], 400.0e6)
        self.assertAlmostEqual(strain[-1], 0.15)
        self.assertGreater(yield_strain, 0.0)


class QuantumTests(unittest.TestCase):
    def test_energy_scales_as_inverse_width_squared(self):
        e1 = energy_levels_ev(1.0e-9, 9.1093837139e-31, 1)[0]
        e2 = energy_levels_ev(2.0e-9, 9.1093837139e-31, 1)[0]
        self.assertAlmostEqual(e1 / e2, 4.0)

    def test_ground_state_wavefunction_is_zero_at_walls(self):
        x = np.array([0.0, 0.5e-9, 1.0e-9])
        psi = stationary_wavefunction(1.0e-9, 1, x)
        self.assertAlmostEqual(psi[0], 0.0, places=12)
        self.assertGreater(psi[1], 0.0)
        self.assertAlmostEqual(psi[2], 0.0, places=12)

    def test_wavefunction_is_normalized(self):
        width = 1.0e-9
        x = np.linspace(0.0, width, 10001)
        psi = stationary_wavefunction(width, 3, x)
        self.assertAlmostEqual(float(np.trapezoid(np.abs(psi) ** 2, x)), 1.0, places=5)


class CrystalTests(unittest.TestCase):
    def test_conventional_cell_basis_counts(self):
        expected = {
            "Simple Cubic (SC)": 1,
            "Body-Centred Cubic (BCC)": 2,
            "Face-Centred Cubic (FCC)": 4,
        }
        for name, count in expected.items():
            with self.subTest(structure=name):
                points = unit_cell_points(name, 3.0e-10)
                self.assertEqual(len(points), count)
                self.assertEqual(STRUCTURES[name]["atoms_per_cell"], count)

    def test_packing_fraction_order(self):
        sc = STRUCTURES["Simple Cubic (SC)"]["packing_fraction"]
        bcc = STRUCTURES["Body-Centred Cubic (BCC)"]["packing_fraction"]
        fcc = STRUCTURES["Face-Centred Cubic (FCC)"]["packing_fraction"]
        self.assertLess(sc, bcc)
        self.assertLess(bcc, fcc)

    def test_density_is_positive(self):
        density = theoretical_density_kg_m3(
            "Face-Centred Cubic (FCC)", 3.61e-10, 63.546e-3
        )
        self.assertGreater(density, 0.0)
        self.assertAlmostEqual(atomic_radius_m("Face-Centred Cubic (FCC)", 3.61e-10), 3.61e-10 / (2.0 * np.sqrt(2.0)))


if __name__ == "__main__":
    unittest.main()
