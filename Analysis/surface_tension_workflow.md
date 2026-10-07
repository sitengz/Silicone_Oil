# Standard film surface-tension analysis (MATLAB)

Use this same procedure for each completed oil film. Keep the original energy
and log files together, and keep different runs of the same case separate.
The analyzer does not submit/stop simulations or modify their files.

## Inputs and calculation

For a case such as `N4_PDI1`, supply:

- `energy.N4_PDI1.film_eq.dat`: film equilibration.
- `energy.N4_PDI1.film.dat`: film production.
- `out.N4_PDI1.film`: matching LAMMPS log (recommended for completion checks).

Both energy clocks start at zero. Production is offset by the final
equilibration time for the trend plot; the shared boundary snapshot is counted
once. Initial film preparation before the reset clock is not part of this axis.
The parser accepts the old nine-column and new eleven-column named headers,
but rejects malformed rows, missing samples/irregular cadence, decreasing or
duplicate times, and mismatched boundary snapshots.

For each saved sample, calculate the two-interface mechanical estimate:

```text
gamma [mN/m] = 0.0101325 * Lz [angstrom] / 2
             * (Pzz - (Pxx + Pyy)/2) [atm]
```

Use the full box `Lz` with the global box-averaged pressure, including vacuum;
do not substitute a density-derived film thickness. If `Lz` varies, average
the per-sample product rather than multiply separately averaged terms. Bulk
energy is not needed for this pressure-tensor method. It is not the alternative
bulk-versus-film energy comparison method.

## Fixed default settings

| Item | Setting |
| --- | --- |
| Trend window | 5 ns |
| Window shift | 1 ns; values plotted at window centers |
| Production uncertainty blocks | Non-overlapping 5 ns |
| Default reported interval | Entire available production stage |
| Late-stage check | Final up-to-20 ns of the selected interval |
| Mean/error bar | Mean of complete production blocks ± block SEM |
| Additional interval | Approximate Student-t 95% CI of block means |
| Figure style | 800 × 600 px window, Times New Roman 14, no title |
| Saved formats | TIFF 300 dpi, PNG preview, editable MATLAB FIG, CSV tables |

All intervals are half-open `[start, end)`. The final endpoint is excluded;
only complete blocks enter the reported mean and SEM. Any incomplete tail is
excluded from block statistics and its duration is listed in the summary. The
all-production and selected-interval sample means are also saved separately.
For uniform sampling and complete blocks, the block mean equals the mean of
the samples used in those blocks.

For `B` complete block means `g_b`, the error bar is
`std(g_b, 0) / sqrt(B)`. It is SEM, not the block SD or a 95% CI. The optional
95% CI is `mean(g_b) ± t_(0.975, B-1) * SEM`; no Statistics Toolbox is needed.
Fewer than two blocks gives no meaningful SEM/CI. A 5 ns block is a starting
setting, not proof that block means are independent.

## One case

Run this in MATLAB after adding the repository's `Analysis` folder to the path:

```matlab
addpath('C:/Users/siteng/Documents/ChatGPT/Si oil descriptor/Silicone_Oil/Analysis');
dataDir = '//wsl.localhost/Ubuntu/home/siteng/OilOnly';
outDir = fullfile(pwd, 'surface_tension_results');
r = analyze_film_surface_tension('N4_PDI1', dataDir, outDir);
disp(r.summary);
```

On Linux MATLAB, set `dataDir = '/home/siteng/OilOnly'` or the actual results
directory. With no `outputDir`, results go to `Analysis/results/<case>/`.
Each case gets full-stage and production-only sliding plots, a production
mean ± SEM figure, and windows/blocks/late-blocks/summary CSV files. Repeating
the same case and output directory replaces its analysis outputs, not inputs.

## Many cases

```matlab
cases = ["N4_PDI1", "Eg1_PDMS_N32_M3125"];
[summary, failures] = analyze_film_surface_tension_batch(cases, dataDir, outDir);
% Or discover all energy.*.film.dat files in this directory:
[summary, failures] = analyze_film_surface_tension_batch([], dataDir, outDir);
```

The batch continues if a case fails validation, recording the reason in
`surface_tension_failures.csv`. Successful cases are collected in
`surface_tension_summary.csv` and a comparison plot with block-SEM error bars.
A successfully analyzed case is NOT automatically a converged or finished run.
Results from earlier cases not included in a later batch remain on disk;
the aggregate CSV represents only the cases processed in that batch.

## Template verification

The included MATLAB tests check nine- and eleven-column inputs, a known
production mean and SEM, interval selection and incomplete tails, invalid
names/missing files/duplicate clocks, and batches containing failed cases:

```matlab
tests = runtests('test_film_surface_tension.m');
assertSuccess(tests);
```

On 2026-10-01 all seven tests passed. The template also reproduced all 96
window values from the original N4 case-specific MATLAB script and analyzed
both N4 and the saved Eg1 data successfully. This verifies the computation,
not convergence of either film.

## Acceptance checks, every time

1. Confirm the log belongs to these energy files, the production stage is
   finished, and there is no `ERROR`. A `Total wall time:` marker is a useful
   completion indicator, not proof of equilibrium. A newer incomplete log
   is flagged because it may belong to a resubmission with older energy files.
2. Inspect wall contact. Nonzero selected-production wall forces require
   review; the calculated value should not be called unconfined-film tension.
   Old files without these columns cannot establish zero wall contact.
3. Inspect both trend plots and non-overlapping blocks. Compare first/last
   half block means, the block slope, and the late-stage mean. Drift is NOT
   included in the statistical SEM, and a small SEM does not establish a plateau.
4. Check SEM sensitivity with 5 and 10 ns blocks (and longer blocks if enough
   data are available); short/serially correlated blocks can understate errors.
5. Confirm the film has two flat, separated interfaces normal to z, does not
   bridge vacuum or hit walls, and the box/pressure conventions are appropriate.
   The energy files alone cannot validate that geometry.
6. If relaxation persists, extend production or report the value as provisional.
   Do not automatically call the whole-stage mean an equilibrium surface tension.

For a later plateau, select the interval explicitly in **production-local** ns
and save it under a different output directory to retain the original analysis:

```matlab
rLate = analyze_film_surface_tension('N4_PDI1', dataDir, ...
    fullfile(outDir, 'late_interval'), 'ProductionStartNs', 30, 'ProductionEndNs', 50);
r10 = analyze_film_surface_tension('N4_PDI1', dataDir, ...
    fullfile(outDir, 'block10'), 'BlockNs', 10);
```

Record the selected interval and block size with every quoted result. The
sliding-window points overlap and are correlated: never use their scatter or
count to calculate the reported SEM.

## N4_PDI1 reference check

The saved N4 data analyzed on 2026-10-01 have 50 ns equilibration followed by
50 ns production, giving 96 overlapping windows. Entire production gives
`55.3038258 ± 0.6980470 mN/m` (SEM, ten complete 5 ns blocks), with approximate
95% CI `[53.7247337, 56.8829179]`. Production-local 30–50 ns gives
`53.9753555 ± 0.2604348 mN/m` (four blocks).

These are provisional: the late trend decreases, and the newer incomplete
log does not establish that the resubmitted job has finished. This example
demonstrates the checks, not a final accepted equilibrium measurement.
