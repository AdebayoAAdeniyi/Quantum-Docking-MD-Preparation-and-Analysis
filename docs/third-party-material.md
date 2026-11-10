# Third-party material: what it is, where to get it, what to do with it

Some of my original working folders held code and parameter files written by
other people. None of it is redistributed here, because it is not mine to
redistribute. This document covers each item: what it was for, which script
touches it, where it comes from, and what to do with it once you have it.

Earlier versions of this repository listed these items without saying where to
get them. That was the gap this document fills.

Links were checked in November 2026.

| item | needed for | status |
|---|---|---|
| [MGLTools](#mgltools-and-autodocktools) | module 13, ligand and receptor PDBQT | replaced by Meeko, not needed |
| [PyReweighting](#pyreweighting) | module 08, aMD reweighting | needed if you run accelerated MD |
| [CHARMM36 for GROMACS](#charmm36-force-field-for-gromacs) | module 09 route B | needed only for that route |
| [Tian Lu's VMD scripts](#tian-lus-vmd-scripts) | module 05, isosurface rendering | optional, for figures |
| [funnel metadynamics tutorial files](#funnel-metadynamics-tutorial-files) | module 12 | optional, reference material |
| [DSSP](#dssp) | nothing | not a dependency; explained below |
| [the R tarballs](#the-r-tarballs) | module 14 | superseded by CRAN and conda |
| [Rosetta](#rosetta) | module 15 | `clean_pdb.py` and `octopus2span.pl` needed for some steps |

## MGLTools and AutoDockTools

**What it was.** My original docking preparation called `prepare_ligand4.py`
and `prepare_receptor4.py` from MGLTools to write the PDBQT files that
AutoDock Vina reads. MGLTools ships its own Python interpreter.

**What this repository does instead.** `13-docking/vina/vina-prepare.sh` calls
Meeko's `mk_prepare_ligand.py` and `mk_prepare_receptor.py`. Meeko is from the
Forli lab and is one of the two preparation routes the Vina documentation
describes, the other being the ADFR suite, which provides `prepare_ligand` and
`prepare_receptor` with its own self-contained Python
([Vina software requirements](https://autodock-vina.readthedocs.io/en/stable/docking_requirements.html)).
Meeko installs with pip or conda, which is why it is the route here.

**If you want MGLTools anyway.** The old Scripps address
`ccsb.scripps.edu/mgltools/` now redirects to `sannerlab.github.io/software/mgltools/`.
Install it in its own environment, not the one holding Meeko, and call the
scripts through the MGLTools Python rather than the system one:

```bash
$MGLROOT/bin/pythonsh $MGLROOT/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py \
    -l ligand.pdbqt -A hydrogens
```

One thing from the original scripts not to carry over: they patched metal
charges by running `sed` on the finished PDBQT. That depends on the atom order
in the file and fails silently when it changes. Metal handling belongs in the
receptor preparation step.

## PyReweighting

**What it is for.** Accelerated MD applies a boost potential, so the resulting
ensemble is not canonical. To get a free energy surface you have to reweight
the trajectory. `08-amber/run/amber-amd-inputs.sh` writes the aMD input with
`-amdlog`, which produces the weight file, and PyReweighting does the
reweighting.

**Where it comes from.** The Miao lab at UNC, repository `MiaoLab20/pyreweighting`
([PyReweighting page](https://www.med.unc.edu/pharm/miaolab/resources/PyReweighting)).
It provides `PyReweighting-1D.py`, `-2D.py` and `-3D.py` with Python 3, bash
wrappers `reweight-1d.sh`, `reweight-2d.sh` and `reweight-3d.sh`, and a
Python 2.7 copy of the 1D and 2D scripts. It implements three reweighting
methods: exponential average, Maclaurin series expansion and cumulant
expansion to second order.

**How it joins this workflow.** After an aMD run, extract the reaction
coordinate with `cpptraj` and the weights from the aMD log, then reweight.
Cumulant expansion is the method the Miao 2014 paper recommends for
reasonable boost magnitudes:

`amber-amd-inputs.sh` writes `run-amd.sh`, which runs the production in chunks
and passes `-amdlog amd-$i.log` to each one, so the boost log you need is
`amd-1.log`, `amd-2.log` and so on.

```bash
# 1. the reaction coordinate, one value per frame, from the aMD trajectory
cpptraj -p solv.prmtop <<'EOF'
trajin amd-1.nc
trajin amd-2.nc
distance D1 :LIG@C1 :315@CA out cv.dat
EOF
awk 'NR>1 {print $2}' cv.dat > Psi.dat

# 2. the weights, from the aMD logs
#    PyReweighting wants one line per frame holding
#      dV/kT   dV   dV_dihedral   dV_total
#    The aMD log carries the dihedral and total boost in separate columns.
#    Check the column positions against the header of your own amd-*.log
#    before using them: the layout has changed between Amber versions, and
#    silently taking the wrong column produces a smooth, wrong free energy
#    surface rather than an error.
head -5 amd-1.log                       # read the header, then build the awk

# 3. reweight, cumulant expansion to second order
python PyReweighting-1D.py -input Psi.dat -T 300 -cutoff 10 -disc 6 \
    -Emax 20 -job amdweight_CE -weight weights.dat
```

The wrappers take positional arguments instead:
`./reweight-1d.sh $Emax $cutoff $binx $data $T`. The `-job` choices are
`amdweight_CE` for cumulant expansion, `amdweight_MC` for the Maclaurin series
and `amdweight` for the exponential average.

Check the boost distribution before trusting the result. If the boost spans
more than about 10 to 20 kT, reweighting is dominated by a few frames and the
free energy surface is noise however it is computed. Gaussian accelerated MD
exists because of this problem.

**Cite** Miao, Sinko, Pierce, Bucher, Walker and McCammon, *J. Chem. Theory
Comput.* **2014**, 10, 2677.

## CHARMM36 force field for GROMACS

**What it was.** A `charmm36-mar2019.ff` directory sat beside my original
GROMACS scripts, which built the protein with CHARMM36 and the ligand with
GAFF2 and then edited the ligand's 1-4 electrostatic scaling so the two could
be combined. That is not a valid combination, and `09-gromacs/README.md`
explains why the module uses the AMBER family on both sides instead.

**Where it comes from.** The MacKerell lab GROMACS ports, at
[CHARMM force field files](https://mackerell.umaryland.edu/charmm_ff.shtml).
The exact file my scripts used, `charmm36-mar2019.ff.tgz`, is still listed
there among the older packages. The current ones at the time of writing are
the February 2026 release, offered as `charmm36-feb2026_cgenff-5.0.ff.tgz` and
`charmm36-feb2026_ljpme_cgenff-5.0.ff.tgz`, with CGenFF 4.6 variants alongside.
The page notes that both CGenFF versions are provided because "the CGenFF
parameter file version must match" the CGenFF program version you use, and
that version 5.0 is the default.

**What to do with it.** Unpack it next to your system so `pdb2gmx` finds it,
then use route B of the module, which is for protein-only systems where
CHARMM36 is used as published:

```bash
tar xzf charmm36-feb2026_cgenff-5.0.ff.tgz    # creates charmm36-*.ff/
gmx pdb2gmx -f protein.pdb -o processed.gro -p topol.top -ff charmm36-feb2026 -water tip3p
```

**If you need CHARMM for a protein with a ligand,** the correct route is CGenFF
for the ligand, not GAFF2. Get the ligand stream file from the CGenFF service,
then convert it with the `cgenff_charmm2gmx_py3_nx2.py` script on the same
MacKerell page. Do not mix CGenFF and GAFF2 either.

## Tian Lu's VMD scripts

**What they are for.** `05-wavefunction-analysis/run/multiwfn-batch.sh` drives
Multiwfn to write cube files for orbitals, electrostatic potential surfaces and
Fukui functions. Turning those cubes into publication figures is done in VMD,
and Multiwfn's author distributes VMD scripts such as `orbitalshow.vmd` for
exactly that.

**Where they come from.** They ship with the Multiwfn manual and example
package from [Tian Lu's site](http://sobereva.com/multiwfn). They are not
redistributed here, and neither are the large saved VMD state files that were
in my original folder.

**How they join this workflow.** Multiwfn writes the cube, VMD renders it:

```bash
# the task names are fukui, mesp, homo-lumo, overlap, hole-electron, custom
05-wavefunction-analysis/run/multiwfn-batch.sh -t homo-lumo -l list   # writes the cubes
vmd -e orbitalshow.vmd                                                # renders one
```

Open the `.vmd` script and set the file name and isosurface value at the top
before running it; the scripts are written to be edited.

**Cite** Lu and Chen, *J. Comput. Chem.* **2012**, 33, 580, for Multiwfn, and
Humphrey, Dalke and Schulten, *J. Mol. Graphics* **1996**, 14, 33, for VMD.

## Funnel metadynamics tutorial files

**What they were.** My original folder held the TCL and shell scripts from the
funnel metadynamics tutorial, set up for one specific system, and they depended
on PLUMED's optional `funnel` module being compiled in.

**Where they come from.** The reference input for the funnel module is on
PLUMED-NEST at [eggs/19/039](https://www.plumed-nest.org/eggs/19/039/), linked
from the [funnel module documentation](https://www.plumed.org/doc-v2.9/user-doc/html/_f_u_n_n_e_l_m_o_d.html).
That page describes `FUNNEL_PS` as the collective variable measuring the
ligand's position relative to a line, and `FUNNEL` as the bias that applies the
restraint, a cone over the binding site opening into a cylinder towards the
solvent. It also points to a full protocol paper by Raniolo and Limongelli
covering the pre- and post-processing stages.

**You do not need them for module 12.** `plumed-input.sh` builds the same
restraint from `PROJECTION_ON_AXIS`, `CUSTOM` and `UPPER_WALLS`, which are in
every PLUMED build, and generates it from residue selections for any system
rather than one. Get the tutorial files to compare against a reference
implementation, or if you prefer to use the real `FUNNEL` action, which needs
`--enable-modules=funnel` at configure time.

**Cite** Limongelli, Bonomi and Parrinello, *PNAS* **2013**, 110, 6358.

## DSSP

**It is not a dependency of anything here, and that is worth stating clearly**,
because a `dssp` tarball sat in my original folder and earlier versions of this
repository listed it as third-party material without saying what used it.

Nothing in this repository calls an external DSSP binary:

* Module 08 assigns secondary structure with `cpptraj`'s own `secstruct`
  command, which is built into cpptraj. It writes files named `dssp.dat` and
  `dssp-summary.dat`, which is a naming habit from the method rather than a
  dependency on the program.
* Module 14's `bio3d-analysis.R` does trajectory PCA, RMSF and dynamic
  cross-correlation. It does not call `bio3d`'s `dssp()` function, which is the
  function that would need the external binary.

If you add a secondary structure step that does need it, DSSP is now maintained
as a rewrite with mmCIF support at [PDB-REDO/dssp](https://github.com/PDB-REDO/dssp),
installs as the `mkdssp` binary, and is packaged for Debian and Ubuntu as
`dssp`. `bio3d::dssp()` takes the path to that binary with its `exefile`
argument.

## The R tarballs

My original folder held source tarballs of R and of several R packages,
because the cluster I was using did not have them. Install them from the
distribution instead:

```bash
micromamba install -c conda-forge r-base r-bio3d r-factominer
# or, inside R
install.packages(c("bio3d", "FactoMineR"))
```

bio3d is from the Grant lab, [thegrantlab.org/bio3d](http://thegrantlab.org/bio3d/),
and FactoMineR from [factominer.free.fr](http://factominer.free.fr).

## Rosetta

**Rosetta is now module 15.** The comparative modelling, docking and
constrained relaxation scripts are in
[`../15-rosetta-modelling/`](../15-rosetta-modelling/), with the corrections
made to them listed in that module's README. An earlier version of this
document recorded the absence of this work as an open gap; that gap is closed.

**Where Rosetta comes from.** RosettaCommons states that "for non-commercial
users (academic, non-profits, government) Rosetta is available for free"
([Rosetta download](https://rosettacommons.org/download/)). Two routes are
offered: release versions as precompiled binaries or source from
`downloads.rosettacommons.org`, and the development version built from source
through GitHub. Commercial licences go through the University of Washington
CoMotion office. Set `ROSETTA_BIN` and `ROSETTA_SUFFIX` in
`config/cluster.local.conf`; the binaries carry a build suffix such as
`relax.default.linuxgccrelease`, which differs between installations.

### Two Rosetta scripts that are not redistributed here

Both ship with Rosetta under `main/tools/` and are other people's work.

**`clean_pdb.py`**, by Phil Bradley, Rhiju Das, Michael Tyka, TJ Brunette and
James Thompson of the Baker Lab, with edits by Steven Combs, Sam Deluca and
Jordan Willis of the Meiler Lab. It extracts one chain from a PDB, strips what
Rosetta cannot read and renumbers from 1. It lives in
`main/tools/protein_tools/scripts/` and needs Python 2.

`rosetta-cm-thread.sh` calls it with `-c` when `ROSETTA_TOOLS` is set, and
otherwise uses the template PDB files as given, which is what you want if you
have already trimmed them in a viewer. Rosetta needs each template as a single
chain numbered continuously from 1, and it will not tell you when that is
untrue; it will thread badly instead.

**`octopus2span.pl`**, by Bjorn Wallner and modified by Vladimir
Yarov-Yarovoy. It converts a predicted transmembrane topology from the OCTOPUS
or SPOCTOPUS server into the `.span` format the Rosetta membrane score terms
read. It lives in Rosetta's membrane tools directory.

`rosetta-cm-setup.sh` calls it when given `-O` and `ROSETTA_TOOLS` is set. Only
membrane proteins need it. Without a `.span` file the membrane weight sets
cannot be used, and modelling a transmembrane protein with the soluble weights
compacts the helices into a globule, because the radius of gyration term is
doing the wrong job; see
[`../15-rosetta-modelling/wts/README.md`](../15-rosetta-modelling/wts/README.md).

### Weight files that are included

The six score function weight files in
[`../15-rosetta-modelling/wts/`](../15-rosetta-modelling/wts/) are included
rather than cited, because a protocol is not reproducible without the score
function it used, and they are term names and numbers rather than program code.
Their provenance, the published soluble RosettaCM weights plus the membrane
adaptation, is set out in that directory's README.

### Web services the module depends on

Three inputs come from servers and cannot be scripted here: fragment libraries
from the Robetta fragment server, transmembrane topology from TOPCONS or
OCTOPUS, and the template alignment from a multiple aligner such as Clustal
Omega. The glycosylation and membrane embedding that followed in the original
work used the GLYCAM glycoprotein builder and CHARMM-GUI, also web services.
[`running.md`](running.md#i-building-a-structure-with-rosetta) records where
each sits in the sequence.
