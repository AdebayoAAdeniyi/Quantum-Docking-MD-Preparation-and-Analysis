# Troubleshooting notes

Problems that actually came up in this work, with what fixed them. They are
collected from the notes kept beside the original scripts.

## Gaussian

| message or symptom | cause and fix |
|---|---|
| `Optimization stopped. -- Number of steps exceeded` | the geometry is still moving; restart from the last structure (`01-gaussian/run/gaussian-status.sh -r`) and raise `MaxCycles`, or optimise in redundant internal coordinates with a better starting structure (xTB pre-optimisation, module 04) |
| `Convergence failure -- run terminated` (SCF) | add `SCF=XQC` (already in the generated inputs) or `SCF=(QC,MaxConventionalCycles=500)`; for metals, start from a smaller basis and read the orbitals with `Guess=Read` |
| `The combination of multiplicity N and M electrons is impossible` | the charge/multiplicity pair is wrong. `gaussian-input.sh` checks this before writing anything |
| `Atoms too close` or `Small interatomic distances` | a bad starting geometry, usually from a 2D structure; rebuild with `00-structure-preparation/build-3d-from-smiles.sh` |
| an imaginary frequency at a minimum | re-optimise with `Opt=(CalcFC)` or `int=ultrafine`; `gaussian-status.sh -r` writes the restart geometries |
| `galloc: could not allocate memory` | `%Mem` exceeds what the node gives; lower it with `-M` |
| NBO or `Pop=NBORead` fails | Gaussian must be linked against NBO; the `nbo` follow-up assumes that |

## ORCA

* **SCF does not converge for a metal complex.** Add `SlowConv` or `VerySlowConv`, use `!NOSOSCF`, or start from the converged orbitals of a smaller basis with `MORead`.
* **CASSCF does not converge, so NEVPT2 is skipped.** Run the first CASSCF with `ExtremeSCF`, then read those orbitals. Check the active space occupations: values close to 0.00 or 2.00 mean the orbital contributes nothing and should be removed. `orca-extract.sh -k casscf` prints a warning for each such orbital.
* **NEB does not converge.** Loosen `Tol_MaxFP_I`, raise `NImages`, or supply a better guess. If the path has to pass a region where atoms clash, build the starting path through a transition state guess rather than by linear interpolation.
* **ESD phosphorescence.** ORCA will not compute the adiabatic S0–T1 gap for you; `orca-esd-input.sh` computes `DELE` from the two optimised states.

## NWChem

* **`no basis on ghost atom`** in a counterpoise job: ghost atoms need their own basis block (`bq<element> library <element> <basis>`), which `nwchem-input.sh -j bsse` writes.
* **NEB starts at an absurd energy.** The straight line between the end points pushed an atom through another one. Supply `NAME-ts` so the starting path is built in two halves.
* **A saddle search walks away from the transition state.** Compute the Hessian first (`task dft freq`) and let the search read it (`inhess 2`), which the generated input does.

## Amber

* **`FATAL: Atom .R<LIG>.A<...> does not have a type`** in tleap: the ligand mol2 and frcmod do not match, or antechamber used a different atom type set. Rebuild both with `amber-ligand-params.sh` (GAFF2 for both).
* **`WARNING: There is a bond of X angstroms between ...`** usually means a missing or wrongly protonated residue in the PDB; check the `pdb4amber` log for gaps.
* **Energy jumps or `vlimit exceeded`** at the start: the minimisation was too short, or the restraint mask is wrong. The two-stage minimisation in `amber-md-inputs.sh` handles most cases.
* **MM-PBSA fails with mismatched atom counts**: the dry topologies must come from the same build as the trajectory. `amber-build.sh` writes complex, receptor and ligand topologies in the same run, so `ante-MMPBSA.py` is not needed.

## GROMACS

* **`Domain decomposition does not support ...` or too few cells.** Use fewer ranks, or `-rdd` / `-dds`; for small systems run with `-ntmpi 1`.
* **`There is no domain decomposition for N ranks`**: the box is too small for the cutoff. Rebuild with a larger `-b`.
* **`Atom X is missing in the topology`** after conversion: the index groups were made from a different structure. Rebuild `index.ndx` from the same `system.gro`.
* **Pressure coupling in 2023**: Berendsen and Parrinello-Rahman were deprecated in favour of C-rescale, which the generated mdp files use.
* **PLUMED**: GROMACS must be patched with PLUMED before `mdrun -plumed` works.

## LAMMPS

* **`Must define pair_style before Pair Coeffs`** when reading a data file written by `write_data`: write it with `write_data ... nocoeff`, or define the styles before `read_data` (the templates include `system.in` first).
* **Restarting**: a data file made from a restart loses the hybrid style coefficients; re-include `system.in` and re-read the coefficients.
* **Green-Kubo viscosity does not plateau**: the run is too short or the correlation length `${p}` is too small; check `viscosity.png` before trusting the number.

## General

* **Open Babel adds hydrogens at the wrong protonation state.** Use `-p 7.4` (`convert-structures.sh -p`) and check charged groups by hand.
* **Different atom order between reactant, product and transition state** breaks QST3, NEB and RMSD comparisons. The scripts refuse to write such inputs.
