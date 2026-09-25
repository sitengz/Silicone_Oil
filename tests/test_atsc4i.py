"""Checks for CG sequence extraction and atomistic ATSC4i reconstruction."""

import importlib.util
import pathlib
import sys
import tempfile
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from Analysis import atsc4i  # noqa: E402


class ATSC4iWorkflowTest(unittest.TestCase):
    def test_capped_smiles_templates(self):
        self.assertEqual(
            atsc4i.capped_smiles("D"), "C[Si](C)(C)O[Si](C)(C)C"
        )
        self.assertEqual(
            atsc4i.capped_smiles("M"), "C[Si](C)(c1ccccc1)O[Si](C)(C)C"
        )
        self.assertEqual(
            atsc4i.capped_smiles("DMD"),
            "C[Si](C)(C)O[Si](C)(c1ccccc1)O[Si](C)(C)O[Si](C)(C)C",
        )
        with self.assertRaises(ValueError):
            atsc4i.capped_smiles("DX")

    def test_extracts_one_cg_chain_in_generator_order(self):
        data = """LAMMPS data file

6 atoms
5 atom types

Atoms # full

1 1 1 0 0 0 0 0 0 0
2 1 4 0 0 0 0 0 0 0
3 1 1 0 0 0 0 0 0 0
4 1 5 0 0 0 0 0 0 0
5 2 1 0 0 0 0 0 0 0
6 2 1 0 0 0 0 0 0 0

Bonds

1 1 1 2
2 1 2 3
3 3 2 4
4 1 5 6

Angles
"""
        with tempfile.TemporaryDirectory() as directory:
            path = pathlib.Path(directory) / "data.test"
            path.write_text(data)
            self.assertEqual(atsc4i.sequence_from_data(path, 1), "DMD")
            self.assertEqual(atsc4i.sequence_from_data(path, 2), "DD")
            with self.assertRaisesRegex(ValueError, "no backbone"):
                atsc4i.sequence_from_data(path, 3)

    def test_config_resolves_case_folder(self):
        with tempfile.TemporaryDirectory() as directory:
            root = pathlib.Path(directory)
            config = root / "model.conf"
            config.write_text("length = 15\noutput = data.07\n")
            path, settings = atsc4i.data_path_from_config(config)
            self.assertEqual(path, root / "data.07")
            self.assertEqual(settings["length"], "15")

    @unittest.skipUnless(
        importlib.util.find_spec("rdkit") and importlib.util.find_spec("mordred"),
        "RDKit and Mordred are not installed",
    )
    def test_mordred_uses_capped_atomistic_graph(self):
        pdms = atsc4i.calculate_atsc4i("D")
        pmps = atsc4i.calculate_atsc4i("M")
        self.assertEqual(pdms["formula"], "C6H18OSi2")
        self.assertEqual(pmps["formula"], "C11H20OSi2")
        self.assertIsInstance(pdms["ATSC4i"], float)
        self.assertIsInstance(pmps["ATSC4i"], float)
        self.assertAlmostEqual(pdms["ATSC4i"], -11.136270936408227, places=6)
        self.assertAlmostEqual(pmps["ATSC4i"], -25.313019676127322, places=6)
        self.assertAlmostEqual(
            atsc4i.calculate_atsc4i("D" * 30)["ATSC4i"],
            -652.3421706239324,
            places=6,
        )


if __name__ == "__main__":
    unittest.main()
