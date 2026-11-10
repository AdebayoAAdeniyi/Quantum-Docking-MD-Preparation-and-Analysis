# Installing the programs these workflows drive

Nothing in this repository needs to be compiled. It is shell, awk, gnuplot and R.
What it needs is the simulation programs themselves, and for some workflows a
few pieces of third-party code that cannot be redistributed here.

Three documents cover this:

* **this file** says how to obtain and configure each program;
* [`../docs/third-party-material.md`](../docs/third-party-material.md) covers the
  code and force field files that are not in the repository, item by item;
* [`../docs/running.md`](../docs/running.md) gives the complete command sequences.

Links were checked in November 2026. Download locations and version numbers
move, so treat the version numbers below as what was current then, not as a
requirement.

## Contents

1. [Free programs, one command](#1-free-programs-one-command)
2. [What each program is used for](#2-what-each-program-is-used-for)
3. [Licensed and registration-required programs](#3-licensed-and-registration-required-programs)
4. [Programs built from source](#4-programs-built-from-source)
5. [Cluster configuration](#5-cluster-configuration)
6. [Checking the installation](#6-checking-the-installation)
7. [What you need for which module](#7-what-you-need-for-which-module)

## 1. Free programs, one command

Everything that was tested while building this repository installs from
conda-forge:

```bash
micromamba create -n qdm -c conda-forge \
    openbabel xtb nwchem gromacs lammps cp2k plumed \
    ambertools acpype parmed \
    vina meeko rdkit gemmi \
    r-base r-bio3d r-factominer gnuplot
micromamba activate qdm
pip install basis_set_exchange
```

`mamba` or `conda` work in place of `micromamba`. If the solver struggles with
the whole set at once, install the MD engines and the chemistry tools in two
separate environments; nothing in this repository requires them to share one.

Two notes on that line:

* **Vina and Meeko.** The AutoDock Vina documentation gives both a pip and a
  conda route, `pip install -U numpy scipy rdkit vina meeko gemmi prody` or
  `conda install -c conda-forge numpy scipy rdkit vina meeko gemmi autogrid`
  ([Vina software requirements](https://autodock-vina.readthedocs.io/en/stable/docking_requirements.html)).
  Module 13 calls `mk_prepare_receptor.py` and `mk_prepare_ligand.py`, which are
  the Meeko commands.
* **PLUMED.** The conda-forge build does not include the optional `funnel`
  module, and it does not need to: module 12 builds the funnel restraint from
  standard PLUMED actions for exactly that reason. See
  [section 4](#4-programs-built-from-source) if you want the real `FUNNEL`
  action, or a PLUMED-patched GROMACS.

## 2. What each program is used for

| program | what it does here | module | script |
|---|---|---|---|
| [Open Babel](https://openbabel.org) | format conversion, 3D generation, protonation at pH | 00, all | `convert-structures.sh`, `build-3d-from-smiles.sh` |
| [xTB](https://github.com/grimme-lab/xtb) | GFN2 screening and pre-optimisation of whole series | 04 | `xtb-run.sh`, `xtb-extract.sh` |
| [NWChem](https://nwchemgit.github.io) | open-source DFT: optimisation, TS, NEB, counterpoise | 04 | `nwchem-input.sh` |
| [AmberTools](https://ambermd.org/AmberTools.php) | GAFF2 parameters, `tleap`, `sander`, `cpptraj`, `MMPBSA.py`, ParmEd | 08, 09, 10 | all of `08-amber/` |
| [GROMACS](https://www.gromacs.org) | MD, free energy landscapes | 09 | `gmx-build.sh`, `gmx-md-inputs.sh` |
| [CP2K](https://www.cp2k.org) | classical and QM/MM MD from an Amber system | 10 | `cp2k-input.sh` |
| [LAMMPS](https://www.lammps.org) | polymer and materials properties | 11 | `11-lammps/` |
| [PLUMED](https://www.plumed.org) | well-tempered and funnel metadynamics | 12 | `plumed-input.sh`, `plumed-fes.sh` |
| [AutoDock Vina](https://vina.scripps.edu) + [Meeko](https://github.com/forlilab/Meeko) | docking and PDBQT preparation | 13 | `13-docking/vina/` |
| [R](https://cran.r-project.org) with [bio3d](http://thegrantlab.org/bio3d/) and [FactoMineR](http://factominer.free.fr) | property statistics, PCA, trajectory PCA, cross-correlation | 14 | `property-stats.R`, `bio3d-analysis.R` |
| [gnuplot](http://www.gnuplot.info) | every plot | all | `plot-md.sh`, `plot-profile.sh`, `plot-neb.sh`, `bin/qdm-plot-spectra` |
| [Basis Set Exchange](https://www.basissetexchange.org) | GAMESS-format basis sets and ECPs | 03 | `gamess-input.sh` |

AmberTools is free. The full Amber package, which adds `pmemd` and
`pmemd.cuda`, is licensed separately from the same site; the scripts in module
08 fall back to `sander` when `pmemd` is absent, so AmberTools alone is enough
to run the protocol on a small system, just slowly.

## 3. Licensed and registration-required programs

These workflows generate the inputs and parse the output. The programs are not
installed by this repository and most cannot be redistributed. Set their
executable names and module load lines in `config/cluster.local.conf`
(see [section 5](#5-cluster-configuration)).

| program | module | how it is obtained | what this repository does with it |
|---|---|---|---|
| **Gaussian** 09/16 | 01, 07 | Commercial licence from [Gaussian Inc.](https://gaussian.com); most universities hold a site licence | writes multi-step input, parses the log for energies, frequencies, NPA/Hirshfeld/CHelpG charges, TD-DFT, NMR, polarisability, EPR |
| **ORCA** 5/6 | 02, 07 | Free for academic and personal use; the vendor states "ORCA is and will remain free for academic and personal use" and distributes it through a registered account and the [ORCA Forum](https://orcaforum.kofo.mpg.de) ([FACCTs](https://faccts.de/orca/)) | writes input, parses ground and excited states, SOC-TD-DFT, CASSCF/NEVPT2, ESD rates, NEB-TS |
| **GAMESS-US** | 03 | Free to academics after a licence agreement from the [Gordon group](https://www.msg.chem.iastate.edu/gamess/) | writes input with mixed basis sets and ECPs, parses NEDA output |
| **Firefly** (PC GAMESS) | 03 | Free for non-commercial use, by registration, from the [Firefly project](http://classic.chem.msu.su/gran/firefly/) | NEDA runs; the `NAO` keyword correction described in `03-gamess-firefly/README.md` applies here |
| **NBO** | 03 | Licensed from the [NBO program](https://nbo.chem.wisc.edu) at Wisconsin; often linked into Gaussian or GAMESS | NEDA and natural population analysis. **Your licence key goes in `config/cluster.local.conf`, which is gitignored.** A key was present in the original `Flyini2` files and has been removed |
| **Multiwfn** | 05 | Free, from [Tian Lu's site](http://sobereva.com/multiwfn) | driven in batch through its menu sequences for Fukui functions, ESP surfaces, hole-electron analysis |
| **AIMAll** | 05 | Commercial, with a free limited version, from [TK Gristmill](http://aim.tkgristmill.com) | QTAIM bond critical points from a Gaussian `.wfx` |
| **POLYRATE** | 06 | Free to academics on request from the [Truhlar group](https://comp.chem.umn.edu/polyrate/) | the `.fu6` output table reader. **Untested**: no `.fu6` file was available |
| **Schrodinger Glide** | 13 | Commercial, [Schrodinger](https://www.schrodinger.com) | best-pose extraction from multi-pose SDF, and its score feeds consensus ranking |
| **MOE** | 13 | Commercial, [Chemical Computing Group](https://www.chemcomp.com) | same as Glide |
| **Rosetta** | 15 | Free for non-commercial users (academic, non-profit, government); release binaries or source from `downloads.rosettacommons.org`, development version on GitHub ([Rosetta download](https://rosettacommons.org/download/)) | comparative modelling, symmetric and local docking, constrained relax. Set `ROSETTA_BIN`, `ROSETTA_SUFFIX`, and `ROSETTA_MPI_SUFFIX` for a parallel build |
| **HADDOCK** | 13 | Free for academic use; the [web server and local version](https://www.bonvinlab.org/software/) are registered through the Bonvin lab | `haddock-air.sh` writes the ambiguous interaction restraint `.tbl` file |

Honest statement of what is tested: the free programs in section 1 were run.
For the programs in this table, the input writers and output parsers were
tested against public reference output, including the cclib test set for
Gaussian 03/09/16, ORCA 4.2/5/6, GAMESS-US 2017 and Firefly 8. Multiwfn and
AIMAll could not be installed, so their drivers are exercised on output
fragments only, and the Multiwfn menu key sequences should be checked once
against your own version before you trust a batch run.

## 4. Programs built from source

Two cases need a build rather than a package.

### PLUMED with the funnel module, and a PLUMED-patched MD engine

The `FUNNEL` and `FUNNEL_PS` actions live in an optional module that, in
PLUMED's words, is "not installed by default" and is enabled by adding
`--enable-modules=funnel` to `./configure`
([PLUMED funnel module](https://www.plumed.org/doc-v2.9/user-doc/html/_f_u_n_n_e_l_m_o_d.html)).

```bash
# PLUMED with the funnel module
./configure --prefix=$HOME/opt/plumed --enable-modules=funnel
make -j8 && make install
source $HOME/opt/plumed/lib/plumed/modulefile   # or set PLUMED_KERNEL

# GROMACS patched for PLUMED, needed for gmx mdrun -plumed
cd gromacs-2024.x
plumed patch -p -e gromacs-2024
# then configure and build GROMACS as usual
```

You do not need either of these for module 12 as written. `plumed-input.sh`
builds the funnel from `PROJECTION_ON_AXIS`, `CUSTOM` and `UPPER_WALLS`, which
are in every PLUMED build, so the conda-forge package is enough. Build the
module if you want the original `FUNNEL` action, and patch GROMACS only if you
drive metadynamics from GROMACS rather than from Amber or CP2K. Amber reads a
PLUMED file through `plumed=1` in `&cntrl` without patching.

### DSSP, if you want the external program

No script here requires it; see
[`../docs/third-party-material.md`](../docs/third-party-material.md). If you
want it anyway, DSSP is now maintained as a rewrite at
[PDB-REDO/dssp](https://github.com/PDB-REDO/dssp), installs the `mkdssp`
binary, and is packaged for Debian and Ubuntu as `dssp`.

## 5. Cluster configuration

Copy the template and edit the copy. The copy is read first and is in
`.gitignore`, so your account details and licence keys never reach a commit:

```bash
cp config/cluster.conf config/cluster.local.conf
$EDITOR config/cluster.local.conf
```

```bash
SCHEDULER=slurm                 # local, pbs or slurm
ACCOUNT=your-allocation
PARTITION=compute
CORES=32
MEM=128G
WALLTIME=24:00:00

GAUSSIAN_SETUP="module load gaussian/g16.c01"
ORCA_SETUP="module load orca/6.0.0; export PATH=\$ORCA_DIR:\$PATH"
AMBER_SETUP="module load amber/22; source \$AMBERHOME/amber.sh"
GROMACS_SETUP="module load gromacs/2024.2; source \$GMXRC"
PLUMED_SETUP="module load plumed/2.9"
NBO_LICENCE_KEY="your-key-here"
```

`SCHEDULER=local` runs jobs in the foreground, which is how the tests run and
how you should try any workflow first.

## 6. Checking the installation

```bash
tests/run-all-tests.sh
```

Each module test checks for the programs it needs and skips with a message
rather than failing when they are absent, so a partial installation gives a
useful report. On a machine with only the section 1 programs, the modules for
Gaussian, ORCA, GAMESS, Multiwfn, AIMAll, Glide, MOE and HADDOCK will skip
their execution checks and still test their parsers.

To test one module:

```bash
tests/test-12-plumed.sh
```

## 7. What you need for which module

If you only want part of this repository, this is the shortest path to each
module.

| you want | install | plus |
|---|---|---|
| structure conversion and 3D building (00) | `openbabel` | |
| fast screening of a series (04) | `xtb` | `openbabel` |
| free-software DFT (04) | `nwchem` | `openbabel` |
| protein-ligand MD end to end (08) | `ambertools` | full Amber for `pmemd.cuda` speed |
| GROMACS MD from an Amber system (09) | `gromacs ambertools parmed` | |
| QM/MM MD (10) | `cp2k ambertools` | |
| materials and polymer properties (11) | `lammps` | |
| metadynamics and binding free energy (12) | `plumed` + an MD engine | see [section 4](#4-programs-built-from-source) |
| docking and consensus ranking (13) | `vina meeko rdkit` | Glide or MOE only for their own scores |
| statistics and trajectory PCA (14) | `r-base r-bio3d r-factominer` | |
| protein structure modelling (15) | nothing from conda | Rosetta, section 3; fragments from the Robetta server |
| Gaussian, ORCA, GAMESS workflows (01, 02, 03) | nothing from conda | the licensed program, section 3 |
