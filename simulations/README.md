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
The original validation examples remain in `examples/` and keep their own
case names under `simulations/`.

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
