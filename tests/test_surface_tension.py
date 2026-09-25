"""Checks for film pressure-anisotropy block averaging."""

import pathlib
import sys
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from Analysis import surface_tension  # noqa: E402


class SurfaceTensionTest(unittest.TestCase):
    def test_config_file_is_adjacent(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = pathlib.Path(directory)
            config = folder / "model.conf"
            config.write_text("output = data.03\n")
            self.assertEqual(
                surface_tension.file_from_config(config, "prod"),
                folder / "energy.03.film.dat",
            )
            self.assertEqual(
                surface_tension.file_from_config(config, "equil"),
                folder / "energy.03.film_eq.dat",
            )

    def test_mechanical_gamma_and_blocks(self):
        row = {"time_fs": 0, "pxx_atm": -20, "pyy_atm": -20,
               "pzz_atm": 0, "lz_A": 300}
        self.assertAlmostEqual(surface_tension.gamma_mn_per_m(row), 30.3975)
        rows = [
            {**row, "time_fs": 0},
            {**row, "time_fs": 5e6, "pxx_atm": -30, "pyy_atm": -30},
            {**row, "time_fs": 10e6, "pxx_atm": -40, "pyy_atm": -40},
        ]
        blocks = surface_tension.block_averages(rows, 5)
        self.assertEqual(len(blocks), 2)
        self.assertEqual(blocks[0][:3], (0, 5, 1))
        self.assertAlmostEqual(blocks[0][3], 30.3975)
        self.assertEqual(blocks[1][:3], (5, 10, 1))
        self.assertAlmostEqual(blocks[1][3], 45.59625)

    def test_rejects_missing_pressure_column(self):
        with tempfile.TemporaryDirectory() as directory:
            path = pathlib.Path(directory) / "energy.film.dat"
            path.write_text("# time_fs pxx_atm pyy_atm lz_A\n0 1 2 3\n5 1 2 3\n")
            with self.assertRaisesRegex(ValueError, "lacks required"):
                surface_tension.read_pressure(path)


if __name__ == "__main__":
    unittest.main()
