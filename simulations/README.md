# Simulations

The numbered folders `01`–`07` contain tracked configurations for the seven
experimental formulations. Run a case from the repository root with, for
example:

```bash
make generate CONFIG=simulations/03/model.conf
```

This writes bulk data, bulk and film inputs, and Slurm scripts into
`simulations/03/`. Generation does not submit either job. The large generated
files and runtime outputs are ignored by Git; each `model.conf` is tracked.
The original validation examples remain in `examples/` and generate beside
their own config files.

| Case | Formulation | Chain length | Chains | Total repeats |
|---|---|---:|---:|---:|
| 01 | PDMS | 30 | 3333 | 99990 |
| 02 | PMPS | 12 | 8333 | 99996 |
| 03 | PDMS/PMPS, random, 5 mol% MPS | 179 | 559 | 100061 |
| 04 | PDMS/PMPS, random, 10 mol% MPS | 42 | 2381 | 100002 |
| 05 | PDMS/PMPS, random, 10 mol% MPS | 44 | 2273 | 100012 |
| 06 | PDMS/PMPS, random, 10 mol% MPS | 65 | 1538 | 99970 |
| 07 | PDMS/PMPS, random, 50 mol% MPS | 15 | 6667 | 100005 |

These configs use the same initial density (0.1 g/cm³), compression target
(0.8 g/cm³), and force-field workflow as the validation examples.

## PDMS chain-length/PDI study

The named `N*_PDI*` folders contain 24 additional pure-PDMS `model.conf`
files: six monodisperse lengths (4, 8, 16, 32, 64, 128) and six non-unit
PDI targets (1.05 through 1.30 in steps of 0.05) at each of Mn = 16, 32,
and 64. Every config gives an explicit `chain_count = length count` list.
The generator has no PDI option or distribution fitting. It only builds the
listed molecules. The N16/PDI1.30 case uses lengths 4–34; other N16 cases
are capped at 32. The N32 and N64 caps are 64 and 128. All cases have at most
100,000 DMS beads and exact integer-repeat Mn.

Run a case from the repository root with, for example:

```bash
make generate CONFIG=simulations/N32_PDI1.2/model.conf
```

This produces initial bulk data, bulk and film LAMMPS inputs, and submit
scripts beside that config. It does **not** submit simulation jobs. Generated
simulation files are ignored by Git; the explicit configs are retained.

The four summary tables are in `table_PDI1.csv`, `table_N16.csv`,
`table_N32.csv`, and `table_N64.csv`, with a readable copy in
[`pdms_series_tables.md`](pdms_series_tables.md). The three 300-dpi TIFF
figures are `SZ_N16.tif`, `SZ_N32.tif`, and `SZ_N64.tif`. Each shows six solid
SZ curves computed from the table's `SZ k` and `SZ rate`, with the actual
integer `model.conf` counts as same-color scatter points. A black dashed line
marks Mn; six case-color dashed lines mark the realized Mw values. Each
histogram is contiguous over its occupied range and unimodal; a small
two-chain adjustment preserves the exact Mn while matching the target PDI.
To regenerate the configs and tables:

```bash
make pdms-series
```

On Linux with librsvg and cairo installed, regenerate the TIFF previews with
`make pdms-figures`.

For the MATLAB versions of the three figures, run:

```matlab
addpath('simulations');
plot_pdms_series
```

The MATLAB script reads the config histograms and table SZ parameters, then
saves editable `.fig` and 300-dpi `.tif` files under
`simulations/figures/`. MATLAB is not required to build the simulation cases.
