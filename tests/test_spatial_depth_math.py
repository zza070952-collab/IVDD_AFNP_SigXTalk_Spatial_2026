"""Pure technical synthetic smoke tests; these are not research results."""
import importlib.util
from pathlib import Path
import unittest
import numpy as np

PATH = Path(__file__).resolve().parents[1] / "scripts/spatial_depth_check.py"
SPEC = importlib.util.spec_from_file_location("depth_check", PATH)
DEPTH = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(DEPTH)


class MoleculeSamplingTest(unittest.TestCase):
    def test_exact_without_replacement_and_stable_spot_seed(self):
        original = np.array([5, 70, 200, 31, 4], dtype=np.int64)
        one = DEPTH.sample_counts(original, 100, DEPTH.stable_rng(20261022, "capture", "spot"))
        two = DEPTH.sample_counts(original, 100, DEPTH.stable_rng(20261022, "capture", "spot"))
        self.assertTrue(np.array_equal(one, two))
        self.assertEqual(int(one.sum()), 100)
        self.assertTrue(np.all(one <= original))
        self.assertTrue(np.all(one >= 0))

    def test_exact_depth_is_unchanged(self):
        original = np.array([10, 20, 70])
        sampled = DEPTH.sample_counts(original, 100, DEPTH.stable_rng(1, "x", "y"))
        self.assertTrue(np.array_equal(sampled, original))

    def test_below_target_is_not_imputed(self):
        with self.assertRaises(ValueError):
            DEPTH.sample_counts(np.array([10, 20]), 100, DEPTH.stable_rng(1, "x", "y"))


if __name__ == "__main__":
    unittest.main()
