"""Check the explicit PDMS chain-count series without generating large data files."""

import csv
import importlib.util
import pathlib
import struct
import unittest
import xml.etree.ElementTree as ET


ROOT = pathlib.Path(__file__).resolve().parents[1]
SIMULATIONS = ROOT / "simulations"
MODULE_SPEC = importlib.util.spec_from_file_location(
    "build_pdms_series", SIMULATIONS / "build_pdms_series.py"
)
SERIES = importlib.util.module_from_spec(MODULE_SPEC)
MODULE_SPEC.loader.exec_module(SERIES)


def read_counts(path):
    rows = {}
    for line in path.read_text().splitlines():
        if line.startswith("chain_count = "):
            _, values = line.split("=", 1)
            length, count = (int(value) for value in values.split())
            if length in rows:
                raise AssertionError(f"Duplicate length in {path}")
            rows[length] = count
    return rows


class PDMSCaseSeriesTest(unittest.TestCase):
    def test_all_24_cases_have_exact_mn_and_valid_pdi(self):
        expected = []
        for n in SERIES.MONO_LENGTHS:
            expected.append((SERIES.case_name(n, 1.0), n, 1.0, n))
        for mean in SERIES.MEANS:
            for pdi in SERIES.TARGETS:
                upper = 34 if mean == 16 and pdi == 1.30 else 2 * mean
                expected.append((SERIES.case_name(mean, pdi), mean, pdi, upper))

        self.assertEqual(len(expected), 24)
        for name, mean, target, upper in expected:
            with self.subTest(case=name):
                config = SIMULATIONS / name / "model.conf"
                self.assertTrue(config.is_file())
                counts = read_counts(config)
                self.assertTrue(counts)
                self.assertGreaterEqual(min(counts), 4)
                self.assertLessEqual(max(counts), upper)
                self.assertTrue(all(count > 0 for count in counts.values()))
                if target > 1.0:
                    self.assertEqual(sorted(counts), list(range(min(counts), max(counts) + 1)))
                    self.assertTrue(SERIES.is_unimodal(list(counts.values())))
                chains, beads, realized_mean, _, realized_pdi = SERIES.statistics(counts)
                expected_chains = ((100_000 if target == 1.0 else 99_999)
                                   // mean)
                self.assertEqual(chains, expected_chains)
                self.assertEqual(beads, mean * chains)
                self.assertLessEqual(beads, 100_000)
                self.assertEqual(realized_mean, mean)
                self.assertLess(abs(realized_pdi - target), 0.0001)
                self.assertIn(f"output = data.{name}", config.read_text())

    def test_four_tables_and_three_tiff_figures_match_cases(self):
        with (SIMULATIONS / "table_PDI1.csv").open() as handle:
            mono = list(csv.DictReader(handle))
        self.assertEqual(len(mono), 6)
        for row in mono:
            counts = read_counts(SIMULATIONS / row["Case"] / "model.conf")
            self.assertEqual(counts, {int(row["N"]): int(row["Chains"])})
            self.assertEqual(sum(n * count for n, count in counts.items()),
                             int(row["Beads"]))
        for mean in SERIES.MEANS:
            with (SIMULATIONS / f"table_N{mean}.csv").open() as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 6)
            self.assertEqual([row["Target PDI"] for row in rows],
                             [f"{target:.2f}" for target in SERIES.TARGETS])
            self.assertTrue(all(float(row["SZ k"]) > 0 for row in rows))
            self.assertTrue(all(float(row["SZ theta"]) > 0 for row in rows))
            tiff = SIMULATIONS / f"SZ_N{mean}.tif"
            with tiff.open("rb") as handle:
                self.assertEqual(handle.read(4), b"II*\x00")
                ifd_offset = struct.unpack("<I", handle.read(4))[0]
                handle.seek(ifd_offset)
                tag_count = struct.unpack("<H", handle.read(2))[0]
                tags = {}
                for _ in range(tag_count):
                    tag, kind, count, value = struct.unpack("<HHII", handle.read(12))
                    tags[tag] = (kind, count, value)
                self.assertEqual(tags[256], (4, 1, 2240))
                self.assertEqual(tags[257], (4, 1, 1680))
                self.assertEqual(tags[259], (3, 1, 1))
                self.assertEqual(tags[262], (3, 1, 2))
                handle.seek(tags[258][2])
                self.assertEqual(struct.unpack("<HHH", handle.read(6)), (8, 8, 8))
            self.assertEqual(tiff.stat().st_size - tags[273][2], tags[279][2])
            svg = ET.parse(SIMULATIONS / f"SZ_N{mean}.svg").getroot()
            ns = {"svg": "http://www.w3.org/2000/svg"}
            self.assertFalse(any("PDMS SZ fits and chain counts" in (node.text or "")
                                 for node in svg.findall("svg:text", ns)))
            # Six solid fitted curves, six scatter series, and Mn + six Mw refs.
            self.assertEqual(len(svg.findall("svg:polyline", ns)), 6)
            self.assertEqual(len(svg.findall("svg:line[@stroke-dasharray]", ns)), 7)
            self.assertLessEqual(len(svg.findall("svg:line[@stroke='#dedede']", ns)), 6)
            labels = [node.text for node in svg.findall("svg:text", ns)]
            self.assertIn("n", labels)
            self.assertIn("M(n)", labels)
            math_labels = [node for node in svg.findall("svg:text", ns)
                           if node.text in ("n", "M(n)")]
            self.assertTrue(all(node.attrib.get("font-style") == "italic"
                                for node in math_labels))
            self.assertFalse(any("Mn (shared)" in (label or "") or
                                 "case Mw" in (label or "") for label in labels))
            expected_points = sum(len(read_counts(SIMULATIONS / row["Case"] / "model.conf"))
                                  for row in rows)
            self.assertEqual(len(svg.findall("svg:circle", ns)), expected_points + 6)


if __name__ == "__main__":
    unittest.main()
