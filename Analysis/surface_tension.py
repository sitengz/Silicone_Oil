#!/usr/bin/env python3
"""Block-average film surface tension from LAMMPS pressure-tensor output."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

ATM_ANGSTROM_TO_MN_PER_M = 0.0101325
REQUIRED_COLUMNS = ("time_fs", "pxx_atm", "pyy_atm", "pzz_atm", "lz_A")


def file_from_config(config: Path, phase: str) -> Path:
    output = None
    for raw in config.read_text().splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        key, separator, value = line.partition("=")
        if separator and key.strip().replace("-", "_") == "output":
            output = Path(value.strip().strip("\"'"))
    if output is None or not output.name:
        raise ValueError(f"{config} needs an output = data.<case> setting")
    data = output if output.is_absolute() else config.parent / output
    case = data.name.removeprefix("data.")
    suffix = ".film_eq.dat" if phase == "equil" else ".film.dat"
    return data.parent / f"energy.{case}{suffix}"


def read_pressure(path: Path) -> list[dict[str, float]]:
    with path.open() as handle:
        header_line = handle.readline().strip()
        if not header_line.startswith("# "):
            raise ValueError(f"{path} has no column header")
        columns = header_line[2:].split()
        if any(column not in columns for column in REQUIRED_COLUMNS):
            raise ValueError(f"{path} lacks required pressure, time, or Lz columns")
        records = []
        for line_number, raw in enumerate(handle, 2):
            if not raw.strip() or raw.lstrip().startswith("#"):
                continue
            fields = raw.split()
            if len(fields) != len(columns):
                raise ValueError(f"{path}:{line_number} has {len(fields)} fields, expected {len(columns)}")
            record = dict(zip(columns, map(float, fields)))
            if any(not math.isfinite(value) for value in record.values()):
                raise ValueError(f"{path}:{line_number} contains a non-finite value")
            if record["lz_A"] <= 0:
                raise ValueError(f"{path}:{line_number} has non-positive Lz")
            if records and record["time_fs"] <= records[-1]["time_fs"]:
                raise ValueError(f"{path}:{line_number} does not advance in time")
            records.append(record)
    if len(records) < 2:
        raise ValueError(f"{path} needs at least two pressure records")
    return records


def gamma_mn_per_m(row: dict[str, float]) -> float:
    lateral = 0.5 * (row["pxx_atm"] + row["pyy_atm"])
    return 0.5 * row["lz_A"] * (row["pzz_atm"] - lateral) * ATM_ANGSTROM_TO_MN_PER_M


def block_averages(rows: list[dict[str, float]], block_ns: float) -> list[tuple[float, float, int, float]]:
    if not math.isfinite(block_ns) or block_ns <= 0:
        raise ValueError("block_ns must be positive and finite")
    width_fs = block_ns * 1.0e6
    # fix print emits both endpoints; omit the final one from half-open bins.
    samples = rows[:-1]
    blocks: dict[int, list[float]] = {}
    for row in samples:
        block = int(math.floor(row["time_fs"] / width_fs))
        blocks.setdefault(block, []).append(gamma_mn_per_m(row))
    return [
        (index * block_ns, (index + 1) * block_ns, len(values), sum(values) / len(values))
        for index, values in sorted(blocks.items())
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--config", type=Path, help="Model config beside generated files")
    source.add_argument("--file", type=Path, help="Film pressure time-series file")
    parser.add_argument("--phase", choices=("equil", "prod"), default="prod")
    parser.add_argument("--block-ns", type=float, default=5.0)
    args = parser.parse_args()

    path = args.file if args.file else file_from_config(args.config, args.phase)
    rows = read_pressure(path)
    blocks = block_averages(rows, args.block_ns)
    print(f"file: {path}")
    print("time_ns_start time_ns_end samples gamma_mN_per_m")
    for start, end, count, mean in blocks:
        print(f"{start:g} {end:g} {count} {mean:.6f}")
    samples = rows[:-1]
    overall = sum(gamma_mn_per_m(row) for row in samples) / len(samples)
    print(f"overall_sample_mean_mN_per_m: {overall:.6f}")
    print("Block values must be inspected for drift; this tool does not certify equilibration.")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as error:
        raise SystemExit(f"Error: {error}") from None
