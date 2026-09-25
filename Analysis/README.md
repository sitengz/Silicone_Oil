# Analysis

`atsc4i.py` calculates **raw** Mordred ATSC4i for one CG oil chain. It reads
the generated initial `data.<case>` file, extracts the type-1 (DMS) and type-4
(MPS) backbone sequence for one molecule, and reconstructs the atomistic
siloxane graph. It does not use CG coordinates or the bulk/film trajectory.

The provisional terminal model is:

```text
Me-[Si(R)(R')-O]n-Si(Me)3
DMS: R = Me, R' = Me
MPS: R = Me, R' = phenyl
```

Thus the first repeat's Si bonds to one terminal methyl; the last repeat's O
bonds to a terminal `Si(Me)3`. Chain order runs from the Si-side end toward
the O-side end, following increasing backbone atom IDs in the generated data.
This convention should be checked against the collaborator's exact drawing
before comparing numerical values.

Install the Python dependencies in your chosen environment (a virtual
environment is recommended):

```bash
python3 -m pip install -r Analysis/requirements.txt
```

From the repository root, after generating a case:

```bash
make atsc4i CONFIG=simulations/03/model.conf MOLECULE=1
```

The tool prints JSON with the selected D/M sequence, capped atomistic SMILES,
formula, raw ATSC4i, and package versions. RDKit may display the canonical
SMILES from the opposite end; the `sequence` field retains generator order.
Pure-chain configs also require the generated data file when used through
`make`. To check a specified D/M
sequence without CG data, run:

```bash
python3 Analysis/atsc4i.py --sequence DMDM
```

`mps_percent` and the random seed in a config do not uniquely specify a
particular random chain's sequence because the generator also uses its RNG
for packing; the generated data file is the source of truth. The generator
numbers backbone beads in sequence order, so a data file with reordered atom
IDs cannot be oriented safely by this tool.
For numerical comparison, also confirm that the collaborator uses the same
atomic-weight table, explicit-hydrogen convention, and descriptor implementation.

## Film surface tension

`surface_tension.py` reads wall-free film pressure output and calculates the
two-surface mechanical estimate from the normal and lateral pressures:

```text
gamma = (Lz / 2) * [Pzz - (Pxx + Pyy) / 2]
```

Run it after a film job to inspect nonoverlapping 5 ns blocks:

```bash
make surface-tension CONFIG=simulations/03/model.conf PHASE=equil
make surface-tension CONFIG=simulations/03/model.conf PHASE=prod
```

`BLOCK_NS=10` changes the block width. The tool prints mN/m and an overall
sample mean. It does not decide when the film is equilibrated; inspect the
block trend and film geometry before treating a value as final. The
`film_eq` and `film` files each start their time counter at zero.
