# Examples

These three configurations keep chain length and repeat count fixed at
`3125 × 32 = 100000` repeat units. The bead count changes because a PMPS
repeat contributes two beads while a PDMS repeat contributes one.

| Configuration | PDMS:PMPS repeats | Beads |
|---|---:|---:|
| `pdms_n32/model.conf` | 100:0 | 100000 |
| `pmps_n32/model.conf` | 0:100 | 200000 |
| `random_50_50_n32/model.conf` | 50:50 | 150000 |

From the repository root, compile the generator and run one config:

```bash
g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic \
    Generator/oil_generator.cpp -o Generator/oil_generator
./Generator/oil_generator --config examples/pdms_n32/model.conf
```

The output goes beside each example's `model.conf`. That directory contains the
initial bulk data, bulk and film inputs, separate Slurm scripts, and a pair
submission script. To queue both stages, enter that directory and run
`bash submit.<case>.pair.sh`. The film job waits for the bulk job's successful
completion and reads its 300 K NPT-equilibrated snapshot.

These are production-scale examples; test smaller overrides such as
`--chains 40 --output data.pilot_pdms` before submitting the full
systems. The root `model.conf` remains a smaller editable example.
