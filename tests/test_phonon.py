import unittest

import numpy as np

from physicsx.phonon import (
    ATOMIC_MASS_UNIT_KG,
    angular_frequency_to_thz,
    diatomic_chain_dispersion,
    phonon_density_of_states,
)


class PhononTests(unittest.TestCase):
    def setUp(self):
        self.result = diatomic_chain_dispersion(
            lattice_constant_m=5.0e-10,
            force_constant_n_m=20.0,
            mass1_kg=23.0 * ATOMIC_MASS_UNIT_KG,
            mass2_kg=35.0 * ATOMIC_MASS_UNIT_KG,
            q_points=801,
        )

    def test_acoustic_mode_is_zero_at_gamma(self):
        q = self.result["q_reduced"]
        index = int(np.argmin(np.abs(q)))
        self.assertAlmostEqual(self.result["acoustic_omega_rad_s"][index], 0.0, places=5)

    def test_optical_mode_is_finite_and_positive_at_gamma(self):
        q = self.result["q_reduced"]
        index = int(np.argmin(np.abs(q)))
        omega = self.result["optical_omega_rad_s"][index]
        self.assertTrue(np.isfinite(omega))
        self.assertGreater(omega, 0.0)

    def test_all_branches_are_real_and_nonnegative(self):
        self.assertTrue(np.all(np.isfinite(self.result["acoustic_omega_rad_s"])))
        self.assertTrue(np.all(np.isfinite(self.result["optical_omega_rad_s"])))
        self.assertTrue(np.all(self.result["acoustic_omega_rad_s"] >= 0.0))
        self.assertTrue(np.all(self.result["optical_omega_rad_s"] >= 0.0))

    def test_equal_masses_touch_at_zone_boundary(self):
        data = diatomic_chain_dispersion(
            lattice_constant_m=4.0e-10,
            force_constant_n_m=12.0,
            mass1_kg=28.0 * ATOMIC_MASS_UNIT_KG,
            mass2_kg=28.0 * ATOMIC_MASS_UNIT_KG,
            q_points=401,
        )
        self.assertAlmostEqual(
            data["acoustic_omega_rad_s"][-1],
            data["optical_omega_rad_s"][-1],
            places=5,
        )

    def test_density_of_states_integrates_to_two_modes_per_cell(self):
        centres, dos, edges = phonon_density_of_states(
            self.result["acoustic_omega_rad_s"],
            self.result["optical_omega_rad_s"],
            bins=120,
        )
        integral = float(np.sum(dos * np.diff(edges)))
        self.assertAlmostEqual(integral, 2.0, places=10)
        self.assertEqual(len(centres), 120)

    def test_frequency_unit_conversion(self):
        converted = angular_frequency_to_thz(np.array([2.0 * np.pi * 1e12]))
        self.assertAlmostEqual(converted[0], 1.0, places=12)


if __name__ == "__main__":
    unittest.main()
