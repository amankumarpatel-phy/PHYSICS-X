import unittest

import numpy as np

from physicsx.crystal_advanced import (
    first_brillouin_zone,
    lattice_geometry,
    miller_d_spacing_m,
    miller_plane_section,
    powder_xrd_pattern,
    reciprocal_lattice_points,
)


class AdvancedCrystallographyTests(unittest.TestCase):
    def test_reciprocal_duality(self):
        for structure in (
            "Simple Cubic (SC)",
            "Body-Centred Cubic (BCC)",
            "Face-Centred Cubic (FCC)",
            "Hexagonal Close-Packed (HCP)",
        ):
            with self.subTest(structure=structure):
                c = 1.633e-10 if structure.startswith("Hexagonal") else None
                geometry = lattice_geometry(structure, 3.0e-10, c)
                direct = geometry["direct"]
                reciprocal = geometry["reciprocal"]
                np.testing.assert_allclose(
                    direct @ reciprocal.T,
                    2.0 * np.pi * np.eye(3),
                    rtol=1e-12,
                    atol=1e-12,
                )

    def test_cubic_d_spacing(self):
        geometry = lattice_geometry("Simple Cubic (SC)", 4.0e-10)
        d = miller_d_spacing_m((1, 1, 1), geometry["reciprocal"])
        self.assertAlmostEqual(d * 1e10, 4.0 / np.sqrt(3.0), places=10)

    def test_miller_plane_section_has_points_on_plane(self):
        geometry = lattice_geometry("Simple Cubic (SC)", 4.0e-10)
        hkl = np.array([1.0, 1.0, 0.0])
        section = miller_plane_section((1, 1, 0), 1.0, geometry["direct"])
        fractional = section @ np.linalg.inv(geometry["direct"])
        np.testing.assert_allclose(fractional @ hkl, 1.0, atol=1e-9)
        self.assertGreaterEqual(len(section), 3)

    def test_bcc_and_fcc_systematic_absences(self):
        bcc = lattice_geometry("Body-Centred Cubic (BCC)", 3.0e-10)
        fcc = lattice_geometry("Face-Centred Cubic (FCC)", 3.0e-10)
        bcc_indices, _ = reciprocal_lattice_points(
            bcc["reciprocal"], "Body-Centred Cubic (BCC)", bcc["basis"], 2
        )
        fcc_indices, _ = reciprocal_lattice_points(
            fcc["reciprocal"], "Face-Centred Cubic (FCC)", fcc["basis"], 2
        )
        self.assertTrue(all((int(h) + int(k) + int(l)) % 2 == 0 for h, k, l in bcc_indices))
        self.assertTrue(all(len({abs(int(v)) % 2 for v in (h, k, l)}) == 1 for h, k, l in fcc_indices))

    def test_brillouin_zone_volume_identity(self):
        for structure in (
            "Simple Cubic (SC)",
            "Body-Centred Cubic (BCC)",
            "Face-Centred Cubic (FCC)",
            "Hexagonal Close-Packed (HCP)",
        ):
            with self.subTest(structure=structure):
                c = 1.633 * 3.0e-10 if structure.startswith("Hexagonal") else None
                geometry = lattice_geometry(structure, 3.0e-10, c)
                bz = first_brillouin_zone(
                    geometry["reciprocal"], structure, geometry["basis"], 3
                )
                expected = (2.0 * np.pi) ** 3 / geometry["primitive_volume_m3"]
                self.assertAlmostEqual(bz["volume_m3"] / expected, 1.0, places=6)

    def test_xrd_has_allowed_reflections_and_finite_curve(self):
        geometry = lattice_geometry("Face-Centred Cubic (FCC)", 3.61e-10)
        x, y, peaks = powder_xrd_pattern(
            geometry["reciprocal"], geometry["basis"], "Face-Centred Cubic (FCC)"
        )
        self.assertTrue(np.all(np.isfinite(x)))
        self.assertTrue(np.all(np.isfinite(y)))
        self.assertGreater(float(np.max(y)), 0.0)
        self.assertGreater(len(peaks), 0)
        # For FCC, allowed reflections have all indices even or all odd.
        for peak in peaks:
            hkl = peak["representative_hkl"].strip("()").split()
            vals = [abs(int(v)) % 2 for v in hkl]
            self.assertEqual(len(set(vals)), 1)


if __name__ == "__main__":
    unittest.main()
