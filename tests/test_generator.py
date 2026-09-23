"""Smoke checks for the bulk-to-film generator package (no LAMMPS needed)."""

import json
import math
import pathlib
import subprocess
import tempfile
import unittest
from collections import Counter


ROOT = pathlib.Path(__file__).resolve().parents[1]
EXAMPLES = (
    ("pdms_n32", 1280, 0),
    ("pmps_n32", 2560, 32),
    ("random_50_50_n32", 1920, 16),
)
FORMULATIONS = (
    ("01", 30, 0, 3333),
    ("02", 12, 100, 8333),
    ("03", 179, 5, 559),
    ("04", 42, 10, 2381),
    ("05", 44, 10, 2273),
    ("06", 65, 10, 1538),
    ("07", 15, 50, 6667),
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

    def test_numbered_formulations_distribute_mps_across_chains(self):
        for case, length, percent, configured_chains in FORMULATIONS:
            with self.subTest(case=case):
                config = ROOT / "simulations" / case / "model.conf"
                settings = dict(
                    line.split("=", 1) for line in config.read_text().splitlines()
                    if "=" in line and not line.lstrip().startswith("#")
                )
                settings = {key.strip(): value.strip() for key, value in settings.items()}
                self.assertEqual(int(settings["length"]), length)
                self.assertEqual(int(settings["chains"]), configured_chains)
                self.assertEqual(float(settings["mps_percent"]), percent)
                self.assertEqual(settings["sequence"], "random")
                self.assertEqual(settings["output"], f"simulations/data.{case}")

                trial_chains = 40
                result = subprocess.run(
                    [
                        str(self.generator), "--config", str(config),
                        "--chains", str(trial_chains), "--output",
                        str(self.work / ("data.form_" + case)),
                    ],
                    cwd=ROOT,
                    text=True,
                    capture_output=True,
                )
                self.assertEqual(result.returncode, 0, result.stderr)
                folder = self.work / ("form_" + case)
                info = json.loads((folder / ("form_" + case + ".info")).read_text())
                composition = info["composition"]
                target_mps = math.floor(length * trial_chains * percent / 100 + 0.5)
                self.assertEqual(composition["total_repeats"], length * trial_chains)
                self.assertEqual(composition["mps_repeats_total"], target_mps)
                self.assertEqual(
                    info["topology_counts"]["atoms"],
                    length * trial_chains + target_mps,
                )
                distribution = composition["chain_composition_distribution"]
                self.assertEqual(sum(part["chains"] for part in distribution), trial_chains)
                self.assertEqual(
                    sum(part["chains"] * part["mps_repeats"] for part in distribution),
                    target_mps,
                )
                self.assertLessEqual(len(distribution), 2)
                data = (folder / ("data.form_" + case)).read_text()
                atoms = data.split("Atoms # full\n\n", 1)[1].split("\nBonds\n", 1)[0]
                actual_mps = Counter()
                for line in atoms.splitlines():
                    if line.strip():
                        fields = line.split()
                        if fields[2] == "5":
                            actual_mps[int(fields[1])] += 1
                actual_counts = Counter(actual_mps.get(i, 0) for i in range(1, trial_chains + 1))
                self.assertEqual(
                    actual_counts,
                    Counter({part["mps_repeats"]: part["chains"] for part in distribution}),
                )

    def test_weight_percent_uses_system_wide_rounding(self):
        result = subprocess.run(
            [
                str(self.generator), "--length", "15", "--chains", "40",
                "--mps-wt", "50", "--output", str(self.work / "data.weight_test"),
            ],
            cwd=ROOT,
            text=True,
            capture_output=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        info = json.loads((self.work / "weight_test/weight_test.info").read_text())
        composition = info["composition"]
        expected = math.floor(600 * 74.0 / (136.2264 + 74.0) + 0.5)
        self.assertEqual(composition["mps_repeats_total"], expected)
        self.assertIsNone(composition["mps_repeats_per_chain"])
        self.assertEqual(
            sum(part["chains"] for part in composition["chain_composition_distribution"]),
            40,
        )


if __name__ == "__main__":
    unittest.main()
