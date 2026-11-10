# 09 · GROMACS: topology, MD and analysis

The same protein-ligand systems as module 08, run in GROMACS. The recommended route converts the Amber system with ParmEd, so the protein and the ligand stay in one force field family and the same system can be run in either program.

```
09-gromacs/
├── prepare/gmx-build.sh      route A: convert an Amber system; route B: pdb2gmx for protein only
├── run/gmx-md-inputs.sh      em, nvt (optional annealing), npt, md + a restartable job script
└── analysis/gmx-analysis.sh  energies, RMSD, RMSF, Rg, H-bonds, SASA, clustering,
                              free energy landscape, protein-ligand interaction energy
```

## Workflow

```bash
# route A (protein + ligand), from module 08
08-amber/prepare/amber-ligand-params.sh -m bcc ligand.sdf
08-amber/prepare/amber-build.sh -n complex -p receptor.pdb -L amber-ligands/ligand
09-gromacs/prepare/gmx-build.sh -n complex -a amber-systems/complex

# route B (protein only)
09-gromacs/prepare/gmx-build.sh -n apo -p protein.pdb -f amber99sb-ildn -w tip3p -b 1.0 -s 0.15

09-gromacs/run/gmx-md-inputs.sh -d gromacs-systems/complex -T 300 -V 200 -N 200 -P 100 -A -g
bin/qdm-submit -s run-md.sh gromacs-systems/complex
09-gromacs/analysis/gmx-analysis.sh -d gromacs-systems/complex -s md -I
```

`gmx-build.sh` writes `index.ndx` with the groups the rest of the workflow expects: `LIG`, `Protein_LIG` and `Water_and_ions`, so the temperature coupling groups and every analysis selection are consistent.


The complete sequence, the same system converted from Amber and run in GROMACS, is in [`../docs/running.md`](../docs/running.md#g-the-same-system-in-gromacs). Every script also prints its own usage block with `-h`.

## The force field decision, and why it differs from the original scripts

The original topology scripts built the protein with **CHARMM36** and the ligand with **GAFF2** through tleap and ParmEd, and then edited the ligand's 1-4 electrostatic scaling from 1.2 to 1.0 so that the two topologies could be combined.

That edit makes the files combine, but the two force fields were parameterised with different rules: different 1-4 scaling, different charge derivation, different water models and different Lennard-Jones combination rules. Mixing them is not a validated combination, and the ligand's intramolecular electrostatics end up scaled differently from the way GAFF2 was fitted.

This module therefore uses the AMBER family on both sides (ff14SB for the protein, GAFF2 for the ligand, the combination both were designed for), which is also what module 08 builds. CHARMM36 is still available through route B for protein-only systems, where it is used as published. If you need CHARMM for the protein, the proper route is CGenFF for the ligand, not GAFF2.

## Analysis

One call produces: potential energy, temperature and density from the `.edr`; backbone and ligand RMSD after `trjconv -pbc mol -center`; per-residue RMSF; radius of gyration; hydrogen bonds; SASA; gromos clustering; and the free energy landscape G = −kT ln(P/P_max) on the RMSD/R_g plane with its contour plot. With `-I` it reruns the trajectory with energy groups to get the protein-ligand Coulomb and Lennard-Jones interaction energies.

Energy terms are selected **by name**, not by the menu number that changes with the system; the original scripts detected those numbers with a chain of greps for each system.

## Corrections made to the original scripts

* The mixed CHARMM36/GAFF2 topology with the hand-edited 1-4 scaling, discussed above.
* **Thirteen topology scripts** (protein only, protein + ligand, glycan, metal, and so on) and **ten run scripts** differing only in cluster modules and trajectory length became one build script and one input generator.
* **Berendsen pressure coupling** was removed in favour of C-rescale, as GROMACS 2021+ recommends; Berendsen does not sample the NPT ensemble.
* **Position restraints.** Systems converted from Amber have no `POSRES` sections, so the restraint stages in the old workflow silently did nothing. The generator now writes a `posre-*.itp` for every solute molecule type and includes it under `#ifdef POSRES`.
* **Restartability**: each stage is skipped if its `.gro` already exists, and `mdrun -cpi` continues an interrupted stage.
* The bundled `charmm36-mar2019.ff` directory is not redistributed. It comes from the MacKerell lab GROMACS ports at [mackerell.umaryland.edu/charmm_ff.shtml](https://mackerell.umaryland.edu/charmm_ff.shtml), where that exact file is still listed among the older packages. [`../docs/third-party-material.md`](../docs/third-party-material.md#charmm36-force-field-for-gromacs) covers how to unpack and use it, and what to do if you need CHARMM for a protein with a ligand.

## Tested

The Abl kinase + imatinib system built in module 08 (43,914 atoms) was converted here with ParmEd, the index groups were generated, position restraints were written for both solute molecule types, and the system ran through energy minimisation and into NVT with GROMACS 2023. `tests/test-09-gromacs.sh` repeats the conversion and the input generation on a small ligand box and runs a short minimisation.
