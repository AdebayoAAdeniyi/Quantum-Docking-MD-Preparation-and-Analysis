# Running the workflows end to end

Each module README shows the few commands for that module. This document shows
the complete sequences, from a starting structure to a number you can put in a
paper, including what to check between stages and what to do when a stage
fails.

Every command below is one you type. Nothing here is pseudocode.

## Contents

* [The pattern every module follows](#the-pattern-every-module-follows)
* [A. Quantum chemistry on a series of molecules](#a-quantum-chemistry-on-a-series-of-molecules)
* [B. Docking with validation](#b-docking-with-validation)
* [C. Protein-ligand MD in Amber](#c-protein-ligand-md-in-amber)
* [D. Binding free energy by MM-GBSA](#d-binding-free-energy-by-mm-gbsa)
* [E. Funnel metadynamics: absolute binding free energy](#e-funnel-metadynamics-absolute-binding-free-energy)
* [F. Accelerated MD and reweighting](#f-accelerated-md-and-reweighting)
* [G. The same system in GROMACS](#g-the-same-system-in-gromacs)
* [H. QM/MM molecular dynamics with CP2K](#h-qmmm-molecular-dynamics-with-cp2k)
* [Jobs, monitoring and restarts](#jobs-monitoring-and-restarts)

## The pattern every module follows

Four steps, in the same order, in every module: **prepare, run, extract, plot.**

```bash
export QDM_ROOT=/path/to/Quantum-Docking-MD-Preparation-and-Analysis
cp config/cluster.conf config/cluster.local.conf      # edit this, it is gitignored
```

Molecules are named in two files rather than inside the scripts:

```bash
cat list                 # one molecule name per line, blank lines and # ignored
1-H
2-OMe
3-NO2

cat chargemult.txt       # name, charge, multiplicity; defaults 0 and 1
1-H    0  1
3-NO2  -1 1
```

Any script run with `-h`, or with no arguments, prints its own usage block.
That block is the reference for its flags; this document shows how the scripts
fit together.

## A. Quantum chemistry on a series of molecules

Screen a series cheaply with xTB, then take the survivors to DFT.

```bash
# 1. build 3D structures from a SMILES table, protonated for pH 7.4
00-structure-preparation/build-3d-from-smiles.sh -i smiles.txt -o prepared \
    -f xyz,pdb,mol2 -p 7.4 -n 30

# 2. GFN2-xTB optimisation of the whole series, then collect the results
04-nwchem-xtb/run/xtb-run.sh -l list -d prepared -o xtb-jobs
04-nwchem-xtb/extract/xtb-extract.sh -l list -p xtb-jobs -o xtb-results

# 3. DFT on the series, multi-step opt+freq then properties
01-gaussian/prepare/gaussian-input.sh -l list -d prepared -o gaussian-jobs \
    -j opt,pop,chelpg -m B3LYP -b "6-311+G(d,p)" -s smd:water
bin/qdm-submit -l list -p gaussian-jobs/

# 4. extract to CSV once the jobs finish
01-gaussian/extract/gaussian-extract.sh -k summary -l list -p gaussian-jobs -o dft
01-gaussian/extract/gaussian-extract.sh -k charges -l list -p gaussian-jobs -o dft

# 5. statistics across the series
Rscript 14-statistics/r-analysis/property-stats.R -i dft/summary.csv -o stats
```

`gaussian-extract.sh` takes `-k KIND`. The kinds are `summary`, `thermo`,
`charges`, `uv`, `freq`, `nmr`, `polar`, `hyperfine`, `jcoupling`, `fukui`,
`geometry`, and `all` for every table at once. Each writes its own CSV into
the `-o` directory, so ask for the ones your question needs.

**What to check.** After step 2, look at the xTB energies before spending DFT
time: a molecule that failed to converge there will usually fail in DFT too.
After step 3, the extractor reports any job without a normal termination and
any optimisation that ended on an imaginary frequency; fix those before using
the CSV, because a transition state hiding in a set of minima will quietly
distort every thermochemistry number downstream.

Replace step 3 with `02-orca/prepare/orca-input.sh` for ORCA, or
`04-nwchem-xtb/prepare/nwchem-input.sh` for NWChem, which is free. The
extractors write the same CSV columns, so steps 4 and 5 do not change.

## B. Docking with validation

Validate the setup by redocking a crystal ligand before docking anything you
care about. This is the cheapest check in computational chemistry and the one
most often skipped.

```bash
# 1. receptor and ligands to PDBQT, box from the crystal ligand plus 6 A
13-docking/vina/vina-prepare.sh -r receptor.pdb -x crystal_ligand.sdf \
    -p 6 -H 7.4 ligands/*.sdf

# 2. dock everything
13-docking/vina/vina-run.sh -e 16 -n 9
#    on a cluster instead:
13-docking/vina/vina-run.sh -e 16 -n 9 -q && bin/qdm-submit -s run-vina.sh .

# 3. scores, best poses as SDF, and redocking RMSD against the crystal pose
13-docking/vina/vina-extract.sh -x crystal_ligand.sdf -o vina-results
```

**The check that matters.** `vina-results/vina-best.csv` carries one row per
ligand, and with `-x` an RMSD column. Include the crystal ligand itself among
the docked ligands. If it does not come back within 2 A of its crystal pose,
stop: the box, the protonation or the receptor preparation is wrong, and every
score in the table is unreliable. Module 13's README records 0.40 A for
imatinib redocked into 1IEP, and 0.90 A for the same molecule rebuilt from
SMILES, which is what a working setup looks like.

Three ways to set the box, depending on what you know:

```bash
-x crystal_ligand.sdf -p 6          # extent of a reference ligand plus padding
-c "12.4 -3.1 26.8" -s "22 22 22"   # explicit centre and size
-a "A:315:CA,A:381:CA"              # midpoint of named atoms
```

Then combine several programs' scores, declaring which ones are better when
high:

```bash
13-docking/consensus/consensus-rank.sh -i all-scores.csv -H "plp:high,gold:high"
```

## C. Protein-ligand MD in Amber

From a docked pose to an analysed trajectory.

```bash
# 1. the ligand: unique atom names and one residue name
00-structure-preparation/ligand-pdb-for-md.sh best-pose.sdf LIG > ligand.pdb

# 2. GAFF2 parameters with AM1-BCC charges
08-amber/prepare/amber-ligand-params.sh -r LIG -o amber-ligands ligand.pdb

# 3. solvated system, plus the dry topologies MM-PBSA needs later
08-amber/prepare/amber-build.sh -n complex -p receptor.pdb \
    -L amber-ligands/ligand -w tip3p -b 10 -B oct -s 0.15

# 4. the six-stage protocol: min1, min2, heat, density, equil, prod
08-amber/run/amber-md-inputs.sh -d amber-systems/complex -T 300 \
    -H 100 -D 100 -E 1000 -P 100 -f 10 -g

# 5. run it
bin/qdm-submit amber-systems/complex

# 6. analysis in one cpptraj pass, with plots
08-amber/analysis/amber-analysis.sh -d amber-systems/complex -y prod.nc -l :LIG -k 5
```

**Check step 2 before going on.** `amber-ligand-params.sh` prints every
`frcmod` line that `parmchk2` had to guess, which it marks `ATTN, need
revision`. Those parameters were invented by analogy. For a routine organic
ligand a handful of such lines is normal; for anything with an unusual
substructure or a metal, they are the parameters to refit before the numbers
mean anything.

**Check step 3.** Read `leap.log`, not just the file list. It reports the net
charge, the ions added and any missing atoms tleap filled in. A wrong net
charge here produces a trajectory that runs perfectly and is wrong.

**Check step 4 before production.** `analysis/energy.dat` from the equilibration
stages should show the density settling near 1.0 g/cm3 and the temperature flat
at the target. A density still drifting at the end of the `density` stage means
production started before the box was equilibrated.

Flags worth knowing: `-g` uses `pmemd.cuda` and is roughly two orders of
magnitude faster than `sander` on a GPU; `-m` sets the restraint mask, which
defaults to everything that is not water or ions; `-c` sets the non-bonded
cutoff.

## D. Binding free energy by MM-GBSA

Run this on the trajectory from walkthrough C. It uses the dry topologies that
`amber-build.sh` already wrote, so `ante-MMPBSA.py` is not needed.

```bash
08-amber/analysis/amber-mmpbsa.sh -d amber-systems/complex -y prod.nc \
    -f 1:last:10 -g 5 -t 15
```

Output: `mmpbsa/summary.csv` with the binding energy and its terms, and
`mmpbsa/top-residues.csv` with the per-residue decomposition. By default the
residues decomposed are those within 6 A of the ligand in the dry complex;
`-r "1-50,100-120"` names them explicitly instead.

Add `-P` for Poisson-Boltzmann as well as generalised Born, which is slower by
a large factor. `-g 8` selects the igb=8 model, which needs mbondi3 radii.

**How to report these numbers honestly.** No entropy term is included, so the
absolute value is not a binding free energy. What MM-GBSA supports is ranking
within a congeneric series run under identical settings. A single compound's
MM-GBSA value compared against an experimental Kd is not a meaningful
comparison, and reviewers increasingly say so.

## E. Funnel metadynamics: absolute binding free energy

This is the workflow that does give an absolute binding free energy. The idea
is a funnel-shaped restraint: a cone over the binding site that opens into a
narrow cylinder reaching into the solvent. The ligand can leave and return,
but the unbound state is confined to a known volume, so the entropy of the
unbound state is defined and the binding constant can be integrated.

It needs an equilibrated system, so start from walkthrough C.

### Step 1: choose the axis and the geometry

The funnel axis runs from a point deep in the pocket, A, to a point out in the
solvent, B. Pick A from residues lining the pocket and B from residues on the
surface beyond the entrance, then read their numbers from the PDB that matches
your MD numbering.

```bash
# look at the pocket: which residues contact the ligand, from walkthrough C
sort -k2 -gr amber-systems/complex/analysis/ligand-contacts.dat | head -20
```

### Step 2: generate the PLUMED input

```bash
ambpdb -p amber-systems/complex/solv.prmtop \
       -c amber-systems/complex/equil.rst7 > complex.pdb

12-plumed-metadynamics/plumed-input.sh -p complex.pdb -l LIG \
    -a "315,381" -b "40,44" \
    -z 1.8 -r 0.1 -A 0.6 -m -0.2 -M 3.0 \
    -s 0.05 -H 1.2 -g 10 -T 300 -P 500 > plumed.dat
```

What those numbers mean, in nanometres and kJ/mol:

| flag | meaning | how to choose it |
|---|---|---|
| `-z 1.8` | where the cone becomes a cylinder | just outside the pocket mouth |
| `-r 0.1` | cylinder radius | small; it sets the unbound volume and enters ΔG directly |
| `-A 0.6` | cone half-angle, radians | wide enough to cover every route out of the site |
| `-m -0.2` `-M 3.0` | walls on the projection | `-M` must reach well into bulk solvent |
| `-s 0.05` | Gaussian width | roughly the resolution you want along the axis |
| `-H 1.2` | initial hill height | a fraction of kT; well-tempering decays it |
| `-g 10` | bias factor | higher explores more and converges more slowly |
| `-P 500` | steps between hills | long enough for the system to respond |

Read the generated `plumed.dat` before running it. It is plain text, and the
atom lists in the `COM` lines are the selections you asked for; if a selection
matched the wrong residues, that is visible here and nowhere later.

### Step 3: run the metadynamics

The same `plumed.dat` drives all three engines.

```bash
# Amber: add two lines to the production &cntrl namelist
sed -i "/&cntrl/a  plumed=1, plumedfile='plumed.dat'," amber-systems/complex/prod.in
cd amber-systems/complex && pmemd.cuda -O -i prod.in -p solv.prmtop \
    -c equil.rst7 -o fmetad.out -r fmetad.rst7 -x fmetad.nc

# GROMACS: needs a PLUMED-patched build (see install/README.md section 4)
gmx mdrun -deffnm md -plumed plumed.dat

# CP2K: the -P flag puts the PLUMED section in the input
10-cp2k/cp2k-input.sh -d amber-systems/complex -j mm -P -n 5000000
```

This is a long run. Funnel metadynamics for a druglike ligand typically needs
hundreds of nanoseconds to converge, and there is no way to shorten that.

### Step 4: watch it while it runs

```bash
# is the ligand actually leaving and coming back?
awk '{print $1, $2}' COLVAR | tail -50

# how fast is the bias still growing?
tail -5 HILLS
```

The projection should range over the whole interval between your walls and
recross it several times. A trace that leaves the site once and stays out, or
never leaves at all, will not give a binding free energy however long it runs.
Hill heights should be decaying; well-tempered metadynamics reduces them as a
region fills.

### Step 5: the profile, the convergence check and ΔG

```bash
12-plumed-metadynamics/plumed-fes.sh -i HILLS -k 2.494339 -n 6 \
    -u 0.5:2.5 -r 0.1 -o fes
```

`-k` is kT in the energy unit of the hills: 2.494339 kJ/mol at 300 K, or
0.596 if your hills are in kcal/mol. `-u SITE_MAX:UNBOUND_MIN` says the bound
site is below 0.5 nm along the projection and the unbound plateau begins at
2.5 nm. `-r` must be the same cylinder radius you used in step 2.

Three outputs:

* `fes.dat`, the free energy profile along the projection;
* `fes-convergence.dat` and `fes.png`, the profile rebuilt from increasing
  fractions of the hills;
* `binding.csv`, with the binding constant and ΔG.

**Read the convergence plot before the number.** The script reconstructs the
profile from 1/6, 2/6 and so on of the hills and plots them together. If the
last two curves still differ in the region that matters, the run is not
converged and `binding.csv` is a number without a meaning. This is the main
failure mode of metadynamics, which is why that plot is produced every time
rather than on request.

The flatness of the unbound plateau is the second check. It is the reference
level that ΔG is measured against, so a plateau that still slopes means the
cylinder has not been sampled evenly.

For a well-tempered run without a funnel, use `-c distance`:

```bash
12-plumed-metadynamics/plumed-input.sh -p complex.pdb -l LIG -a "315,381" \
    -c distance -M 2.5 -s 0.05 -H 1.2 -g 10 -T 300 -P 500 > plumed.dat
```

That gives a profile along a distance but no absolute ΔG, because the unbound
volume is not defined.

## F. Accelerated MD and reweighting

aMD adds a boost potential to flatten barriers, so it samples more but no
longer samples the canonical ensemble. Recovering a free energy surface needs
reweighting, which is where PyReweighting comes in; see
[`third-party-material.md`](third-party-material.md) for where to get it.

```bash
# 1. boost parameters computed from a finished conventional run
08-amber/run/amber-amd-inputs.sh -d amber-systems/complex -o prod.out -n 500 -T 300

# 2. run it; run-amd.sh does 100 ns chunks and passes -amdlog to each
bin/qdm-submit -s run-amd.sh amber-systems/complex
```

The thresholds follow the Amber aMD tutorial: `EthreshP = <EPtot> + 0.16 N_atoms`
with `alphaP = 0.16 N_atoms`, and `EthreshD = <DIHED> + 4 N_res` with
`alphaD = 0.8 N_res`, the averages taken from the conventional run's output.
This is why step 1 needs a finished cMD run and will refuse without one.

Then reweight, as set out in
[`third-party-material.md`](third-party-material.md#pyreweighting).

**The honest caveat.** Check the spread of the boost before trusting any
reweighted surface. If it spans more than roughly 10 to 20 kT, the exponential
average is dominated by a few frames and the result is noise regardless of
which reweighting method you choose. Gaussian accelerated MD was developed
because of this, and it is the better choice for new work.

## G. The same system in GROMACS

Route A converts the Amber system from walkthrough C, so the protein and ligand
keep one consistent force field family and the same system can run in both
engines.

```bash
# 1. convert the tleap system with ParmEd
09-gromacs/prepare/gmx-build.sh -n complex -a amber-systems/complex -o gmx-systems/complex

# 2. em, nvt, npt, md
09-gromacs/run/gmx-md-inputs.sh -d gmx-systems/complex -T 300 -V 100 -N 100 -P 100 -f 10 -g
bin/qdm-submit gmx-systems/complex

# 3. analysis, including the free energy landscape and interaction energy
09-gromacs/analysis/gmx-analysis.sh -d gmx-systems/complex -s md -b 1000 -k 0.1 -I
```

Route B is for a protein with no ligand, built with `pdb2gmx`:

```bash
09-gromacs/prepare/gmx-build.sh -n apo -p protein.pdb -f amber99sb-ildn -w tip3p -b 1.0
```

**Two things the original scripts got wrong here**, both covered in
`09-gromacs/README.md`. Position restraints silently did nothing on topologies
converted from Amber, because those carry no `POSRES` sections; `gmx-build.sh`
now generates them with `gmx genrestr` and includes them under `#ifdef POSRES`.
And the protein was built with CHARMM36 while the ligand used GAFF2, combined
by hand-editing the 1-4 scaling, which is not a valid combination. If you need
CHARMM for the protein, use CGenFF for the ligand; see
[`third-party-material.md`](third-party-material.md#charmm36-force-field-for-gromacs).

## H. QM/MM molecular dynamics with CP2K

Treat the ligand, or a reacting subsystem, with DFT and the rest with the Amber
force field, starting from the same system as walkthrough C.

```bash
# classical first, to confirm the conversion is sound
10-cp2k/cp2k-input.sh -d amber-systems/complex -j mm -n 10000 -T 300 -E NPT_I
bin/qdm-submit -s run-cp2k.sh amber-systems/complex

# then QM/MM: the ligand in DFT, electrostatic embedding
10-cp2k/cp2k-input.sh -d amber-systems/complex -j qmmm -q LIG -Q 0 \
    -n 10000 -t 0.5 -T 300 -E NVT
```

`-q` takes an Amber residue name or an explicit atom list such as
`"1-25,40"`. `-Q` is the charge of the QM region and must be right. The time
step drops to 0.5 fs for QM/MM, which the script sets by default.

**Run the classical job first.** It is cheap and it is the only way to tell a
broken topology conversion from a QM problem. Two things that used to break
here, both fixed and explained in `10-cp2k/README.md`: force field atom types
are not element symbols, so CP2K needs explicit `&KIND` sections with
`ELEMENT`, and Amber now writes NetCDF restarts, which need converting to PDB
with `cpptraj` before CP2K can read the coordinates.

## Jobs, monitoring and restarts

Every `run` script writes a `run.sh` in the job directory. `qdm-submit` sends
it to whichever scheduler `config/cluster.local.conf` names:

```bash
bin/qdm-submit amber-systems/complex              # one directory
bin/qdm-submit -l list -p gaussian-jobs/          # every NAME in list
bin/qdm-submit -s run-amd.sh amber-systems/complex   # a differently named script
```

With `SCHEDULER=local` the job runs in the foreground. Try any new workflow
that way first, on a small system, before queueing a hundred of them.

**Restarts.** The multi-stage scripts are written to be resumable. Each stage
is guarded by a test on its own output, so re-submitting after a wall-clock
timeout or a node failure skips the stages that finished and resumes at the
first incomplete one. There is nothing to edit and no stage to comment out.

**When a stage fails,** look at the stage's own output file rather than the
scheduler log: `leap.log` for tleap, the `.out` file for Amber,
`docked/NAME.log` for Vina, `receptor.log` for Meeko. The extractors also
report jobs that did not terminate normally rather than quietly leaving them
out of the CSV. `docs/troubleshooting.md` lists the failures that came up most
often and what each one actually means.
