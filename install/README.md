# Installing the programs these workflows drive

Nothing in this repository needs to be compiled: it is shell, awk, gnuplot and R.
What it needs is the simulation programs themselves.

## Free programs (everything here was tested against these)

Most of them install in one line with conda/mamba from conda-forge:

```bash
micromamba create -n qdm -c conda-forge \
    openbabel xtb nwchem gromacs lammps cp2k plumed \
    ambertools acpype parmed \
    autodock-vina meeko rdkit \
    r-base r-bio3d r-factominer gnuplot
```

| program | what it is used for here | module |
|---|---|---|
| Open Babel | any structure format, 3D generation, protonation | all |
| xTB | GFN2 screening, pre-optimisation | 04 |
| NWChem | open-source DFT: optimisation, TS, NEB, counterpoise | 04 |
| AmberTools | GAFF2 parameters, tleap, sander, cpptraj, MMPBSA.py, ParmEd | 08, 09 |
| GROMACS | MD | 09 |
| CP2K | MM and QM/MM MD from Amber topologies | 10 |
| LAMMPS | polymer and materials properties | 11 |
| PLUMED | metadynamics, funnel metadynamics | 12 |
| AutoDock Vina + Meeko | docking | 13 |
| R (bio3d, FactoMineR) | PCA, cross-correlation, statistics | 14 |
| gnuplot | every plot | all |
| Basis Set Exchange (`pip install basis_set_exchange`) | GAMESS basis sets and ECPs | 03 |

## Licensed programs (the workflows generate their inputs and read their output)

Gaussian, ORCA (free for academics, registration required), GAMESS-US, Firefly,
Multiwfn, AIMAll, NBO, Schrodinger Glide, MOE, POLYRATE, HADDOCK. Set the
executable names and the module load lines in `config/cluster.conf`.

## Cluster settings

Copy `config/cluster.conf` to `config/cluster.local.conf` and edit it; that file
is read first and is in `.gitignore`. It holds the scheduler (local, PBS or
SLURM), the account and queue, the resources, and one setup line per program:

```bash
SCHEDULER=slurm
ACCOUNT=def-yourgroup
CORES=32
MEM=128G
GAUSSIAN_SETUP="module load gaussian/g16.c01"
AMBER_SETUP="module load amber/22; source \$AMBERHOME/amber.sh"
```

## Checking the installation

```bash
tests/run-all-tests.sh        # runs every module test that its programs allow
```
