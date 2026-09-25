# Silicone Oil

This repository contains a standalone generator for neutral PDMS oil, neutral
PMPS oil, and coarse-grained PDMS-PMPS copolymer oils. It was split from the
[`Oil` directory of Silicone_Coating](https://github.com/sitengz/Silicone_Coating/tree/main/Oil).
The generator and its configuration parser are self-contained. `Analysis/`
contains a single-chain ATSC4i tool; surface-tension analyzers are not
included yet.

## Repository layout

```text
Generator/    oil generator and configuration parser
Analysis/     single-chain ATSC4i tool and documentation
simulations/  numbered formulation configs and generated cases
examples/     validation example configurations
model.conf    current editable generator configuration
```

The three example configs are in `examples/`. Each generator run writes a
case-named folder under the requested output directory, following the
elastomer generator's layout:

```text
simulations/<case>/
├── data.<case>                 initial bulk data
├── in.<case>                   bulk LAMMPS input
├── in.<case>.film              dependent film LAMMPS input
├── submit.<case>.sh            bulk Slurm job
├── submit.<case>.film.sh       film Slurm job
├── submit.<case>.pair.sh       submit both with afterok dependency
└── <case>.info                 model and workflow metadata
```

Only bulk data exist at generation time. The film job reads the equilibrated
bulk snapshot written during the bulk run, so it inherits the exact molecules,
sequence, and topology rather than building a second random packing.

The `.info` file is valid JSON and records the composition, realized MPS
monomer and weight percentages, type populations, topology counts, box size,
random seeds, mixing-rule choice, and both simulation stages. During the runs,
LAMMPS also writes separate bulk and film energy/pressure time series.

## Build and generate

From the repository root:

```bash
make
make generate CONFIG=examples/pdms_n32/model.conf
```

`make` only compiles the generator; `make generate` creates the configured
bulk-and-film case. The root `model.conf` is used when `CONFIG` is omitted.
For a small trial run, use:

```bash
make generate CONFIG=examples/pdms_n32/model.conf GENERATOR_ARGS="--chains 40 --output simulations/data.pilot_pdms"
```

Run `make test` for the generator smoke tests. `make clean` removes only the
compiled generator, not generated simulation cases. The generator requires
C++17 for case-folder creation. If `make` is unavailable, the equivalent
compile command is:

```bash
g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic \
    Generator/oil_generator.cpp -o Generator/oil_generator
```

To calculate raw ATSC4i for one generated CG chain, see
[`Analysis/README.md`](Analysis/README.md). After installing its Python
dependencies, run `make atsc4i CONFIG=simulations/03/model.conf MOLECULE=1`.

## Configuration file

Edit [`model.conf`](model.conf) and run:

```bash
./Generator/oil_generator --config model.conf
```

The file uses `key = value` lines. Blank lines and `#` comments are ignored;
keys may use underscores or hyphens (for example, `mps_percent` or
`mps-percent`). All command-line options listed below are available as keys,
including `output`. Command-line values override the same keys in the file:

```bash
./Generator/oil_generator --config model.conf --chains 100
```

The included root configuration describes 500 chains of 32 repeat units with
50% MPS. The config path and `output` path are resolved relative to the
directory where the generator is run. Command-line-only usage remains
available. Generating files does not submit jobs or run LAMMPS.

## Default model

Running without options:

```bash
./Generator/oil_generator
```

generates a bulk case containing 625 PMPS oil chains with 16 MPS repeat units
per chain, plus a film input that waits for the equilibrated bulk data:

```text
Oil_PMPS_N16_M625/
├── data.Oil_PMPS_N16_M625
├── in.Oil_PMPS_N16_M625
├── in.Oil_PMPS_N16_M625.film
├── submit.Oil_PMPS_N16_M625.sh
├── submit.Oil_PMPS_N16_M625.film.sh
├── submit.Oil_PMPS_N16_M625.pair.sh
└── Oil_PMPS_N16_M625.info
```

The bulk contains 20,000 beads and reproduces the counts of the
supplied PMPS N16 model:

| Quantity | Count |
|---|---:|
| Atoms | 20,000 |
| Bonds | 19,375 |
| Angles | 27,500 |
| Dihedrals | 35,000 |

Unlike the old PMPS generator, the initial chains are non-collinear, molecules
are randomly rotated, every complete chain is placed wholly inside the primary
box with zero image flags, and intermolecular placements below 4.5 Å are
rejected. The default initial density is deliberately low at 0.1 g/cm³.

## Atom types

The combined type namespace is reserved now so that the oil can later be added
to V22 or V35 without renumbering:

| Atom type | Meaning | Mass (g/mol) | Used by oil-only generator |
|---:|---|---:|---|
| 1 | Neutral DMS bead | 74.0 | Yes |
| 2 | Reactive DMS type from formulation model | 74.0 | Reserved |
| 3 | Reactive DMS type from formulation model | 74.0 | Reserved |
| 4 | MPS backbone bead | 59.1204 | Yes |
| 5 | MPS phenyl pendant bead | 77.106 | Yes |

One DMS repeat contains one type-1 bead and has a mass of 74.0 g/mol.

One MPS repeat contains one type-4 backbone bead plus one type-5 pendant bead:

```text
M(MPS repeat) = 59.1204 + 77.106 = 136.2264 g/mol
```

Types 2 and 3 have zero atoms in an oil-only system. They are retained in the
data-file header and pair matrix for compatibility with the formulation models.

## Temperature-dependent nonbonded interactions

All pair coefficients are calculated by the generator and written explicitly.
LAMMPS does not apply an automatic mixing rule.

The DMS-MPS combining rule is:

```text
epsilon(DMS,MPS) = 0.579966
                   * sqrt[epsilon(DMS,DMS) * epsilon(MPS,MPS)]
sigma(DMS,MPS)   = [sigma(DMS,DMS)   + sigma(MPS,MPS)]   / 2
```

The epsilon prefactor is applied to DMS interactions with both the MPS
backbone and MPS pendant beads, at both 300 K and 800 K.

The pure PMPS backbone-pendant term remains the supplied PMPS cross
interaction; it is not recalculated by the DMS-MPS combining rule.

### 300 K attractive matrix

At 300 K, the input uses:

```lammps
pair_style lj/gromacs 12 15
special_bonds lj 0 0 0.5
```

| Pair | Epsilon (kcal/mol) | Sigma (Å) | Source |
|---|---:|---:|---|
| DMS-DMS | 1.012878450 | 6.445843660 | V22/V35 DMS model |
| DMS-MPS backbone | 0.778707979 | 5.960508330 | 0.579966 × geometric ε; arithmetic σ |
| DMS-MPS pendant | 0.537772752 | 6.206674705 | 0.579966 × geometric ε; arithmetic σ |
| MPS backbone-backbone | 1.779864125 | 5.475173000 | PMPS model |
| MPS backbone-pendant | 1.331357250 | 5.727290650 | Supplied PMPS cross term |
| MPS pendant-pendant | 0.848858275 | 5.967505750 | PMPS model |

Types 1, 2, and 3 receive the same DMS nonbonded parameters. Therefore all
three types have the same combined interactions with types 4 and 5.

### 800 K repulsive matrix

The DMS parameters at 800 K are calculated with the same temperature equations
used by the V22/V35 generators. Because only 300 K PMPS parameters were
supplied, this standalone model applies the corresponding DMS temperature
ratios to each PMPS epsilon and sigma:

```text
epsilon_PMPS(800) = epsilon_PMPS(300)
                    * epsilon_DMS(800)/epsilon_DMS(300)

sigma_PMPS(800)   = sigma_PMPS(300)
                    * sigma_DMS(800)/sigma_DMS(300)
```

After scaling the pure PMPS terms, the mixed DMS-MPS values are recalculated
using `0.579966 ×` geometric epsilon and arithmetic sigma.

| Pair | Epsilon at 800 K | Sigma at 800 K (Å) | Repulsive cutoff (Å) |
|---|---:|---:|---:|
| DMS-DMS | 0.531850637 | 6.640519401 | 7.453731009 |
| DMS-MPS backbone | 0.408890459 | 6.140526096 | 6.892507499 |
| DMS-MPS pendant | 0.282378187 | 6.394127126 | 7.177165031 |
| MPS backbone-backbone | 0.934585852 | 5.640532791 | 6.331283990 |
| MPS backbone-pendant | 0.699080133 | 5.900264835 | 6.622823352 |
| MPS pendant-pendant | 0.445725560 | 6.147734850 | 6.900599052 |

The high-temperature stage uses `pair_style lj/cut`. Every pair receives its
own cutoff at `2^(1/6)*sigma`, so only the repulsive-force side of that pair's
Lennard-Jones interaction is retained. The global cutoff is the largest
pair-specific value, 7.453731009 Å.

Both complete 800 K and 300 K matrices are also recorded in each `.info` file.

## Rough bonded model

The copolymer rule is implemented by inspecting every bonded interaction.

### Bonds

| Bond type | Interaction | Potential |
|---:|---|---|
| 1 | DMS-DMS or mixed DMS-MPS backbone bond | DMS harmonic bond |
| 2 | MPS-MPS backbone bond | PMPS backbone harmonic bond |
| 3 | MPS backbone-pendant bond | PMPS pendant harmonic bond |

### Angles

| Angle type | Interaction | Potential |
|---:|---|---|
| 1 | Backbone-only angle containing at least one DMS repeat | DMS harmonic angle |
| 2 | Three MPS backbone beads | PMPS backbone quartic angle |
| 3 | Any angle containing an MPS pendant | PMPS pendant quartic angle |

### Dihedrals

| Dihedral type | Interaction | Potential |
|---:|---|---|
| 1 | Backbone-only dihedral containing at least one DMS repeat | DMS `nharmonic` |
| 2 | Four MPS backbone beads | PMPS backbone `nharmonic` |
| 3 | Dihedral containing one MPS pendant | PMPS one-pendant `nharmonic` |
| 4 | Dihedral containing two MPS pendants | PMPS two-pendant `nharmonic` |

Thus a mixed bonded interaction without a pendant follows the DMS model. Once
a pendant participates, it follows the corresponding PMPS rule.

The input uses:

```lammps
bond_style harmonic
angle_style hybrid harmonic quartic
dihedral_style nharmonic
```

## Composition controls

The generator rounds the requested system-wide MPS count to the nearest whole
repeat, then distributes it across chains. Every chain has either the lower
or upper neighboring MPS count; which chains receive the extra repeat is
shuffled reproducibly using `seed`. For `random` sequences, MPS positions
within each chain are shuffled independently. The `.info` file records the
actual system-wide composition and the number of chains at each composition.

### Monomer percentage

Use `--mps-percent` to specify the percentage of repeat-unit positions that are
MPS:

```bash
./Generator/oil_generator \
    --length 32 \
    --chains 500 \
    --mps-percent 25 \
    --sequence random
```

This requests 8 MPS and 24 DMS repeat units per chain. When the requested
percentage does not divide evenly across chains, the generator mixes the two
neighboring per-chain counts to match the overall percentage within one
repeat unit.

### Weight percentage

Use `--mps-wt` to choose the closest overall MPS count based on the different
DMS and MPS repeat masses:

```bash
./Generator/oil_generator \
    --length 32 \
    --chains 500 \
    --mps-wt 40 \
    --sequence random
```

For a chain with `N` repeat positions and `k` MPS repeats:

```text
M_chain = (N-k)*74.0 + k*136.2264

MPS wt% = 100 * k*136.2264 / M_chain
```

Because the total number of MPS repeats must be an integer, the realized
weight percentage can differ slightly from the request. The generator reports
both values and records them in `.info`.

`--mps-percent` and `--mps-wt` cannot be supplied together.

## Sequence modes

| Mode | Behavior |
|---|---|
| `random` | Selects each chain's allocated number of MPS sites randomly and reproducibly |
| `alternating` | Distributes the MPS sites as evenly as possible along each chain |
| `block` | Places one contiguous MPS block in the middle of each chain |

Pure PDMS and pure PMPS chains are unaffected by the sequence choice.

Examples:

```bash
# Pure PDMS oil
./Generator/oil_generator --length 32 --chains 500 --mps-percent 0

# Pure PMPS oil
./Generator/oil_generator --length 16 --chains 625 --mps-percent 100

# Approximately alternating 50:50 copolymer
./Generator/oil_generator \
    --length 32 \
    --chains 500 \
    --mps-percent 50 \
    --sequence alternating

# Central PMPS block at a target MPS weight percentage
./Generator/oil_generator \
    --length 64 \
    --chains 250 \
    --mps-wt 30 \
    --sequence block
```

## All command-line inputs

| Option | Type | Default | Meaning |
|---|---|---:|---|
| `--length N` | positive integer | 16 | Repeat-unit positions per chain |
| `--chains M` | positive integer | 625 | Number of oil chains |
| `--n N` | positive integer | 16 | Alias for `--length` |
| `--m M` | positive integer | 625 | Alias for `--chains` |
| `--mps-percent X` | 0–100 | 100 | MPS monomer percentage |
| `--mps-wt X` | 0–100 | unset | Target MPS repeat-unit weight percentage |
| `--sequence MODE` | text | `random` | `random`, `alternating`, or `block` |
| `--density X` | positive number | 0.1 | Initial mass density in g/cm³ |
| `--target-density X` | positive number | 0.8 | Density after scripted 800 K compression |
| `--film-padding X` | positive number | 300 K repulsive cutoff | Minimum vacuum padding added to each z face, in Å |
| `--min-separation X` | positive number below 15 | 4.5 | Minimum intermolecular bead distance in Å |
| `--seed N` | positive integer | 20260727 | Sequence, conformation, rotation, and packing seed |
| `--velocity-seed N` | positive integer | 492845 | LAMMPS initial-velocity seed |
| `--output FILE` | path | automatic | Name the initial data file; files go in a case-named subdirectory |
| `--config FILE` | path | unset | Read settings from a `key = value` file |
| `--help` | — | — | Print command help |

## Bulk-to-film workflow

The bulk input follows the V22/V35 elastomer equilibration sequence, without
crosslinking:

```text
bulk:  1M steps at 800 K → 1M isotropic compression → 1M relaxation
       → 2M more at 800 K → 1M cooling under isotropic NPT
       → 1M at 300 K under isotropic NPT → write data.<case>.npt_eq
       → 20M-step, 100 ns bulk NVT production
film:  read data.<case>.npt_eq → expose two z surfaces
       → 100k steps at 300 K with temporary walls and lateral NPT
       → remove walls → 1M-step wall-free NVT relaxation
       → 5M-step, 25 ns wall-free NVT production
```

The film starts at **300 K** from the equilibrated bulk snapshot; it does not
repeat hot compression. The input keeps the bulk snapshot's `Lx` and `Ly` at
conversion, while lateral NPT during temporary-wall relaxation may change
their final values. It converts z to nonperiodic (`p p f`) and initially adds
one 300 K repulsive-wall cutoff of space to each z face. The default DMS wall
cutoff is about 7.235 Å. If unwrapped chains would otherwise contact the
walls, the script automatically increases the padding and prints the actual
value in the LAMMPS log. Thus:

```text
film cell Lz = equilibrated bulk Lz + 2 × actual padding
```

This is the cell height, not the oil slab thickness. The film conversion
follows LAMMPS's [bulk-to-slab image-flag procedure](https://docs.lammps.org/Howto_bulk2slab.html)
to preserve bonded chains that cross the original periodic z seam. The
repulsive `wall/lj126` fixes are removed **before** the film's free-surface
relaxation and energy measurement. Check the film trajectory for atoms
approaching the fixed z boundaries after wall removal; increase
`film_padding` if needed.

The high-temperature repulsive pair matrix is active only in bulk's first
five million steps. Both 300 K production stages use the same attractive
`lj/gromacs` matrix, 5 fs timestep, and NVT ensemble. The film's shorter
production samples its surface energy and pressure; it does not collect
Green–Kubo stress data for viscosity. Check time-block averages for drift
before treating 25 ns as sufficient sampling.

### Green-Kubo stress output

After the final 300 K NPT equilibration, the barostat is removed and the final
equilibrated volume is held fixed for a 100 ns NVT trajectory at 300 K. The
production timestep counter and accumulated time are reset to zero.

The atom-coordinate dump is stopped before this long production stage. Instead,
LAMMPS writes the instantaneous off-diagonal pressure components every 10
timesteps, corresponding to a 50 fs sampling interval:

```text
gk_stress.<case>.dat
```

The file has one header followed by exactly four columns:

```text
# time_fs pxy_atm pxz_atm pyz_atm
```

Thus a completed 100 ns run contains 2,000,000 stress samples. The pressure
components are in atmospheres because the input uses LAMMPS `real` units. No
time averaging is applied before writing, preserving the instantaneous stress
series needed for Green-Kubo autocorrelation analysis.

Both production stages also write `energy.<case>.bulk.dat` and
`energy.<case>.film.dat` every 1000 steps. Columns are time, temperature,
potential energy, diagonal pressures, and all three box lengths. The film
file is recorded after walls are removed. A difference in mean potential
energy divided by twice the film cross-sectional area is an energy-based
surface excess, **not automatically the thermodynamic surface tension** at
finite temperature. The recorded pressure components allow a separate
  mechanical-route analysis later. No surface analyzer is included yet.

Each generated Slurm job requests 48 hours. Confirm with shorter tests that
the bulk and film runs each fit this limit; the 100k–200k-bead examples may
need longer allocations or shorter pilot runs.

The scripts use the same Nova module configuration as the V22/V35 generators.
From the generated case folder, submit the ordered pair with:

```bash
bash submit.<case>.pair.sh
```

This submits bulk first, then queues film with Slurm `afterok` dependency.
The film submit script checks for the equilibrated bulk data file and exits
if it is absent. To run the stages manually, submit `submit.<case>.sh`, wait
for successful completion, then submit `submit.<case>.film.sh`.

## Scope of this version

This version generates bulk and film inputs with a true bulk-to-film handoff.
It does not:

- generate a substrate-supported film;
- create reactive oil end groups;
- calculate surface tension automatically;
- infer any mixing rule other than the explicitly documented
  `0.579966 ×` geometric-epsilon, arithmetic-sigma DMS-MPS rule.

## Generated file descriptions

- `data.<case>` is the initial LAMMPS data file containing the box,
  oil atoms, and bonded topology.
- `in.<case>` is the bulk equilibration and production input.
- `in.<case>.film` reads the bulk equilibrated data at runtime and performs
  film conversion, temporary-wall initiation, wall removal, and film production.
- `submit.<case>.sh` and `submit.<case>.film.sh` are separate one-node, 96-task
  Slurm jobs; `submit.<case>.pair.sh` queues them in order.
- `<case>.info` is a JSON manifest containing composition, sequence,
  force-field, topology, random-seed, and production settings.
- `gk_stress.<case>.dat` is the bulk stress series for viscosity analysis.
- `energy.<case>.bulk.dat` and `energy.<case>.film.dat` are production energy,
  pressure, and box-size time series for future surface analysis.
