import unittest

import numpy as np

from physicsx.optics import (
    LorentzOscillator,
    drude_lorentz_dielectric,
    kramers_kronig_eps1_finite_window,
    normal_incidence_slab_rt,
    optical_conductivity_s_m,
    optical_constants,
)


class OpticalModelTests(unittest.TestCase):
    def setUp(self):
        self.energy = np.linspace(0.2, 8.0, 1201)
        self.eps_inf = 1.8
        self.oscillators = [
            LorentzOscillator(2.5, 0.35, 0.8),
            LorentzOscillator(4.2, 0.55, 0.35),
        ]
        self.epsilon = drude_lorentz_dielectric(
            self.energy, self.eps_inf, 8.5, 0.4, self.oscillators
        )

    def test_passive_dielectric_has_nonnegative_imaginary_part(self):
        self.assertTrue(np.all(np.imag(self.epsilon) >= 0.0))
        self.assertTrue(np.all(np.isfinite(self.epsilon)))

    def test_optical_constants_are_passive_and_finite(self):
        result = optical_constants(self.energy, self.epsilon)
        self.assertTrue(np.all(result["n"] >= 0.0))
        self.assertTrue(np.all(result["k"] >= 0.0))
        self.assertTrue(np.all(result["absorption_per_m"] >= 0.0))
        self.assertTrue(np.all(np.isfinite(result["wavelength_m"])))

    def test_conductivity_real_part_is_nonnegative_for_passive_model(self):
        sigma = optical_conductivity_s_m(self.energy, self.epsilon, self.eps_inf)
        self.assertTrue(np.all(np.real(sigma) >= -1e-10))
        self.assertTrue(np.all(np.isfinite(sigma)))

    def test_coherent_slab_energy_balance(self):
        optical = optical_constants(self.energy, self.epsilon)
        # Reconstruct the selected passive complex refractive index.
        index = optical["n"] + 1j * optical["k"]
        result = normal_incidence_slab_rt(
            self.energy, index, 320e-9, 1.0, 1.5
        )
        total = result["reflectance"] + result["transmittance"] + result["absorptance"]
        np.testing.assert_allclose(total, np.ones_like(total), atol=1e-10)
        self.assertTrue(np.all(result["reflectance"] >= 0.0))
        self.assertTrue(np.all(result["transmittance"] >= 0.0))
        self.assertTrue(np.all(result["absorptance"] >= 0.0))

    def test_zero_thickness_air_to_air_is_transparent(self):
        index = np.ones_like(self.energy, dtype=complex)
        result = normal_incidence_slab_rt(
            self.energy, index, 0.0, 1.0, 1.0
        )
        np.testing.assert_allclose(result["reflectance"], 0.0, atol=1e-14)
        np.testing.assert_allclose(result["transmittance"], 1.0, atol=1e-14)
        np.testing.assert_allclose(result["absorptance"], 0.0, atol=1e-14)

    def test_lorentz_resonance_increases_dielectric_loss(self):
        single = drude_lorentz_dielectric(
            self.energy,
            eps_inf=1.0,
            plasma_energy_ev=5.0,
            drude_damping_ev=0.5,
            oscillators=[LorentzOscillator(3.0, 0.15, 1.0)],
        )
        resonance_index = int(np.argmin(np.abs(self.energy - 3.0)))
        self.assertGreater(np.imag(single[resonance_index]), 0.0)

    def test_kramers_kronig_estimate_is_finite_interior(self):
        reconstructed = kramers_kronig_eps1_finite_window(
            self.energy, np.imag(self.epsilon), self.eps_inf
        )
        self.assertTrue(np.all(np.isfinite(reconstructed[1:-1])))
        self.assertTrue(np.isnan(reconstructed[0]))
        self.assertTrue(np.isnan(reconstructed[-1]))

    def test_energy_grid_validation(self):
        with self.assertRaises(ValueError):
            drude_lorentz_dielectric(
                np.array([1.0, 0.5] + [1.0 + i for i in range(8)]),
                1.0, 5.0, 0.5, [],
            )


if __name__ == "__main__":
    unittest.main()
