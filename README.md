# Quantum, Docking and Molecular Dynamics: preparation, running and analysis

A working Linux toolkit for computational chemistry: input preparation, job submission, result extraction, analysis and plotting for the quantum chemistry, molecular dynamics and docking programs used in real research projects, plus the processing of the experimental electrochemistry that goes with them.

These are my own scripts, reorganised and corrected. They were used for the calculations behind **more than 70 peer-reviewed research articles** [Github scientific-research-and-publications](https://github.com/AdebayoAAdeniyi/scientific-research-and-publications):
[Google Scholar](http://scholar.google.co.uk/citations?user=QDAoQt8AAAAJ) · [ResearcherID C-8843-2014](http://www.researcherid.com/rid/C-8843-2014) · [Scopus](https://www.scopus.com/authid/detail.uri?authorId=55388336000) · [ORCID 0000-0003-2421-1382](http://orcid.org/0000-0003-2421-1382)

The original scripts did their job but were personal: cluster paths and molecule names were written inside them, each cluster had its own copy, and a transpose program was pasted into about forty of them. This repository keeps what they knew about the science and makes it usable by anyone: one shared library, one cluster configuration file, lists of molecules instead of hardcoded names, CSV output instead of a pile of temporary files, and a test for every module.

---

## What is here

| | module | programs | what it does |
|---|---|---|---|
| 00 | [structure-preparation](00-structure-preparation) | Open Babel | any input format (pdb, xyz, mol2, sdf, SMILES) to any other, 3D generation, conformer search, protonation, splitting multi-structure files |
| 01 | [gaussian](01-gaussian) | Gaussian 09/16 | multi-step jobs (opt+freq, NPA/Hirshfeld, CHelpG, TD-DFT, NMR, polarisability, EPR, wfx), TS/QST3/IRC, counterpoise, thermochemistry at many temperatures, conceptual DFT |
| 02 | [orca](02-orca) | ORCA 5/6 | ground and excited states, TD-DFT with spin-orbit coupling, CASSCF/NEVPT2, fluorescence and phosphorescence rates (ESD), NEB-TS |
| 03 | [gamess-firefly](03-gamess-firefly) | GAMESS-US, Firefly, NBO | natural energy decomposition analysis (NEDA), mixed basis sets and ECPs from the Basis Set Exchange |
| 04 | [nwchem-xtb](04-nwchem-xtb) | NWChem, xTB | free-software DFT: optimisation, transition states, NEB, counterpoise; GFN2-xTB screening of whole series |
| 05 | [wavefunction-analysis](05-wavefunction-analysis) | Multiwfn, AIMAll | Fukui densities, electrostatic potential surfaces, hole-electron analysis, QTAIM bond critical points |
| 06 | [kinetics](06-kinetics) | any of the above | barriers, Eyring rates with Wigner tunnelling, Arrhenius fits, energy profiles, POLYRATE tables |
| 07 | [electrochemistry](07-electrochemistry) | Gaussian/ORCA/NWChem + lab data | reduction potentials and pKa from a Born-Haber cycle; cyclic voltammetry and UV-Vis processing, HOMO/LUMO energies |
| 08 | [amber](08-amber) | AmberTools, Amber | GAFF2 parameters, solvated systems, the full MD protocol, accelerated MD, cpptraj analysis, MM-GBSA with decomposition |
| 09 | [gromacs](09-gromacs) | GROMACS | topologies (converted from Amber or built with pdb2gmx), MD protocol, full analysis, free energy landscape |
| 10 | [cp2k](10-cp2k) | CP2K | classical and QM/MM molecular dynamics starting from an Amber system |
| 11 | [lammps](11-lammps) | LAMMPS | polymer and materials properties: glass transition, elastic modulus, viscosity, structure and diffusion |
| 12 | [plumed-metadynamics](12-plumed-metadynamics) | PLUMED | well-tempered and funnel metadynamics, free energy profiles, binding free energies |
| 13 | [docking](13-docking) | AutoDock Vina, HADDOCK, Glide, MOE | preparation, docking, redocking validation, interaction restraints, consensus ranking |
| 14 | [statistics](14-statistics) | R (FactoMineR, bio3d) | correlations, PCA and regression of molecular properties; trajectory PCA and cross-correlation |
| 15 | [Protein structure modelling with Rosetta](15-rosetta-modelling/) | comparative modelling from several templates, symmetric and local docking, relax under experimental distance restraints |

Shared parts: [`lib/`](lib) (common shell functions and the awk output parsers), [`bin/`](bin) (`qdm-submit`, `qdm-plot-spectra`), [`config/cluster.conf`](config/cluster.conf), [`examples/`](examples) and [`tests/`](tests).

Documentation:

| document | what it covers |
|---|---|
| [`install/README.md`](install/README.md) | how to obtain and configure every program, free and licensed, and what you need for which module |
| [`docs/running.md`](docs/running.md) | the complete command sequences, end to end, with what to check between stages |
| [`docs/third-party-material.md`](docs/third-party-material.md) | the code and force field files that are not in this repository, and where each comes from |
| [`docs/troubleshooting.md`](docs/troubleshooting.md) | the failures that came up most often and what each one means |

Every script also prints its own usage block when run with `-h` or with no arguments.

## How it works

Everything follows the same four steps, with the same conventions:

```bash
# 1. prepare: a list of molecules, not names inside the script
cat list                       #   1-H
                               #   2-OMe
                               #   MoOCl4
cat chargemult.txt             #   MoOCl4 -1 2       (only what is not a neutral singlet)
01-gaussian/prepare/gaussian-input.sh -l list -d structures -j opt,pop,td,nmr -s smd:acetonitrile

# 2. run: the scheduler comes from config/cluster.conf (local, PBS or SLURM)
bin/qdm-submit -l list -p gaussian-jobs/

# 3. check and extract to CSV
01-gaussian/run/gaussian-status.sh -l list -r
01-gaussian/extract/gaussian-extract.sh -k all -l list

# 4. plot
bin/qdm-plot-spectra -t ir -i gaussian-results/freq.csv -s 0.967 -e experiment.dat
```

The same shape holds in every module, so moving from Gaussian to ORCA or from Amber to GROMACS does not mean learning a new set of habits. Outputs are CSV, so they go straight into R, pandas or a spreadsheet, and the modules feed each other: xTB geometries into DFT, DFT thermochemistry into kinetics and electrochemistry, Amber systems into GROMACS and CP2K, docking poses into MD.

### One configuration file for any machine

```bash
cp config/cluster.conf config/cluster.local.conf   # edited locally, never committed
SCHEDULER=slurm; ACCOUNT=def-yourgroup; CORES=32; MEM=128G
GAUSSIAN_SETUP="module load gaussian/g16.c01"
```

The same commands then run on a laptop, a PBS cluster or a SLURM cluster. The original scripts had one copy per machine.

## Verified, not just written

Programs that are free were **run here**, on real molecules, and the results are in `examples/`:

| what was run | result |
|---|---|
| **xTB** screening of five 2-hydroxybenzophenones | geometries, gaps, IP, EA and reactivity indices; all five are minima |
| **NWChem** B3LYP/6-31G\* | water, HCN and HNC optimised; the HCN → HNC transition state found (one imaginary mode at −1147 cm⁻¹); water dimer counterpoise, BSSE 0.87 kcal/mol |
| **Kinetics** of HCN → HNC | ΔG‡ = 45.4 kcal/mol, Wigner κ = 2.28, Arrhenius fit over 250–400 K (R² = 1.000); the computed HNC−HCN energy is 1.2 kcal/mol above the high-level reference, the expected error for this level |
| **AutoDock Vina** on PDB 1IEP (Abl kinase + imatinib) | the crystal ligand redocks to **0.40 Å** RMSD, and the same molecule rebuilt from SMILES to 0.90 Å |
| **AmberTools** | Abl kinase + imatinib built and solvated (43,914 atoms, no tleap errors); a ligand-in-water box run through the whole protocol, temperature and density converging correctly |
| **GROMACS** | the Amber system converted with ParmEd, index groups and restraints generated, energy minimisation and NVT run |
| **CP2K** | classical MD and QM/MM MD (PBE-D3, GEEP embedding) on the solvated Amber system |
| **LAMMPS** | a polymer melt built and run: glass transition T\*g = 0.569, Young's modulus, Green-Kubo viscosity |
| **PLUMED** | well-tempered metadynamics on a double well whose barrier is exactly 5: recovered 5.27 |
| **R** | PCA and regression on the xTB property table |

Programs that are licensed (Gaussian, ORCA, GAMESS, Firefly, Multiwfn, AIMAll, Glide, MOE) could not be run here, so their **parsers are tested against public reference outputs** from the [cclib](https://github.com/cclib/cclib) project: Gaussian 03/09/16, ORCA 4.2/5.0/6.0, GAMESS-US 2017 and Firefly 8. Where no public output exists (ORCA NEB and ESD, NBO's NEDA table, AIMAll critical points), the parsers are tested on short fragments that copy the layout printed in the program manuals, and [`tests/fixtures/README.md`](tests/fixtures/README.md) says exactly which is which.

```bash
tests/run-all-tests.sh          # every module, skipping what is not installed
tests/run-all-tests.sh 08 12    # only these modules
tests/run-all-tests.sh -s       # include the slow QM and MD runs
```

A test that needs a program you have not installed reports `skipped` and does not fail the run, so a partial installation still gives a useful report. Three ORCA and NBO fixtures are currently absent from the repository, excluded at commit time by a `*.out` rule in `.gitignore`; the tests that need them now say so and skip, and [`tests/fixtures/README.md`](tests/fixtures/README.md) explains how to restore them.

## Errors found in the original scripts

Reviewing old scripts honestly means reporting what was wrong with them. These change results, not just style:

* **Conceptual DFT**: the ionisation potential and electron affinity were computed with the wrong sign, and the electrophilicity index as μ²/2 × η instead of μ²/(2η).
* **NMR extraction**: `grep C` also matched chlorine, and `grep H` matched any line containing an H.
* **Multi-step Gaussian jobs**: taking the last match in the log mixed gas phase and solvent values from different steps. Logs are now split by step and every row says which step it came from.
* **Open-shell NPA charges**: NBO prints total, α and β summaries, and the script took the last one, which is β only.
* **NBO fragment lists** in Firefly inputs were written without the `NAO` keyword, so NBO read the atom numbers as molecular units.
* **Mixed force fields**: protein with CHARMM36 and ligand with GAFF2, combined by editing the 1-4 scaling. Replaced by a consistent AMBER-family setup.
* **Consensus docking ranking**: every score was assumed to be better when lower, which silently inverted the ranking of fitness-style scores.
* **GROMACS position restraints** did nothing on topologies converted from Amber, because those have no `POSRES` sections.
* **Berendsen barostats** (Amber and GROMACS) do not sample the correct NPT ensemble; replaced by Monte Carlo and C-rescale.
* **CP2K** stopped on current systems because force field atom types are not element symbols, and because Amber now writes binary restarts.
* **LAMMPS** `write_data` wrote pair coefficients into the data file, so the next run refused to start.

Each module README lists its own corrections with the reasoning.

## Privacy and third-party material

Cluster accounts, project codes, personal paths, email addresses and an NBO licence key that appeared in the original scripts have been removed. The NBO key now belongs in `config/cluster.local.conf`, which is gitignored.

Third-party code that used to sit beside the scripts is **not** redistributed, because it is not mine to redistribute. [`docs/third-party-material.md`](docs/third-party-material.md) covers each item: what it was for, which script touches it, where to download it and what to do with it once you have it. That list is MGLTools, PyReweighting, the CHARMM36 force field directory for GROMACS, Tian Lu's VMD scripts, the funnel metadynamics tutorial files, DSSP and the R tarballs.

One honest gap in the same document: Rosetta demo files were in my original material, but **no script here uses Rosetta**, so this repository currently under-represents that part of the work rather than merely declining to ship someone else's files. Adding it needs the original protocol scripts.

## About me

I hold bachelor's, master's and PhD degrees in **chemistry**, specialising in **computational chemistry**, and I was a university professor for about a decade from 2016. I completed a degree in **Computing Science at the University of Alberta** (June 2026). The scripts here are the computational side of that research career: the quantum chemistry, molecular dynamics, docking and data analysis behind more than 70 published articles.

I now build **FINWISE**, the main service of my company **Strategic Financial Building** ([www.strategicfinancialbuilding.ca](https://www.strategicfinancialbuilding.ca)), applying the same approach (careful method, verified results, reproducible analysis) to financial planning.

## Licence

MIT, for the scripts in this repository. The programs they drive have their own licences, and published work should cite the programs and the methods, not this repository alone.
