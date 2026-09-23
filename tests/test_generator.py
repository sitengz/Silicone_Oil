"""Smoke checks for the bulk-to-film generator package (no LAMMPS needed)."""

import json
import pathlib
import subprocess
import tempfile
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[1]
EXAMPLES = (
    ("pdms_n32", 1280, 0),
    ("pmps_n32", 2560, 32),
    ("random_50_50_n32", 1920, 16),
)


class GeneratorWorkflowTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(prefix="silicone-oil-test-")
        cls.work = pathlib.Path(cls.temp.name)
        cls.generator = cls.work / "oil_generator"
        subprocess.run(
            [
                "g++", "-std=c++17", "-O2", "-Wall", "-Wextra", "-Wpedantic",
                str(ROOT / "Generator/oil_generator.cpp"), "-o", str(cls.generator),
            ],
            check=True,
        )

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def test_three_example_configs(self):
        for example, expected_atoms, mps_per_chain in EXAMPLES:
            with self.subTest(example=example):
                case = example.replace("_n32", "")
                result = subprocess.run(
                    [
                        str(self.generator), "--config",
                        str(ROOT / "examples" / example / "model.conf"),
                        "--chains", "40", "--output",
                        str(self.work / ("data." + case)),
                    ],
                    cwd=ROOT,
                    text=True,
                    capture_output=True,
                )
                self.assertEqual(result.returncode, 0, result.stderr)
                folder = self.work / case
                info = json.loads((folder / (case + ".info")).read_text())
                self.assertEqual(info["topology_counts"]["atoms"], expected_atoms)
                self.assertEqual(
                    info["composition"]["mps_repeats_per_chain"], mps_per_chain
                )
                self.assertEqual(info["composition"]["chain_length"], 32)
                self.assertEqual(info["composition"]["chain_count"], 40)
                self.assertTrue((folder / ("data." + case)).is_file())
                bulk = (folder / ("in." + case)).read_text()
                film = (folder / ("in." + case + ".film")).read_text()
                self.assertIn("boundary        p p p", bulk)
                self.assertIn("write_data      data." + case + ".npt_eq", bulk)
                self.assertIn("read_data       data." + case + ".npt_eq", film)
                self.assertIn("reset_atoms     image all", film)
                self.assertIn("boundary p p f", film)
                self.assertIn("unfix           zlo_wall", film)
                self.assertIn("unfix           zhi_wall", film)
                self.assertLess(
                    film.index("unfix           zhi_wall"),
                    film.index("file energy." + case + ".film.dat"),
                )
                self.assertNotIn("fix             xlink", bulk + film)
                for suffix in (".sh", ".film.sh", ".pair.sh"):
                    subprocess.run(
                        ["bash", "-n", str(folder / ("submit." + case + suffix))],
                        check=True,
                    )

    def test_rejects_nonpositive_padding(self):
        result = subprocess.run(
            [str(self.generator), "--film-padding", "0"],
            cwd=ROOT,
            text=True,
            capture_output=True,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("--film-padding must be positive", result.stderr)


if __name__ == "__main__":
    unittest.main()
