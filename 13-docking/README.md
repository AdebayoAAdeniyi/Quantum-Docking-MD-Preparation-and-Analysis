# 13 · Docking: Vina, HADDOCK, Glide/MOE and consensus ranking

Preparation, docking, pose extraction and redocking validation with AutoDock Vina; ambiguous interaction restraints for HADDOCK; best-pose extraction from the multi-pose output of any commercial program; and the consensus ranking that combines several scoring functions.

```
13-docking/
├── vina/
│   ├── vina-prepare.sh   receptor and ligand PDBQT (Meeko), search box from a reference ligand,
│   │                     explicit centre, or the midpoint of named atoms
│   ├── vina-run.sh       dock every prepared ligand, locally or as a job
│   └── vina-extract.sh   scores, best poses as SDF, redocking RMSD against a crystal pose
├── consensus/consensus-rank.sh   AASS, Z-score, rank-by-rank and exponential consensus ranking
├── haddock/haddock-air.sh        ambiguous interaction restraints (.tbl) from active/passive residues
└── glide-moe/docking-best-poses.sh  best pose per ligand from any multi-pose SDF with a score field
```

## Workflow

```bash
13-docking/vina/vina-prepare.sh -r receptor.pdb -x crystal_ligand.sdf -p 6 -H 7.4 ligands/*.sdf
13-docking/vina/vina-run.sh -e 16
13-docking/vina/vina-extract.sh -x crystal_ligand.sdf      # scores + redocking RMSD
13-docking/consensus/consensus-rank.sh -i all-scores.csv -H "plp:high"
```

The box can be defined three ways: from a reference ligand (`-x`, its extent plus padding), from an explicit centre (`-c "x y z" -s "22 22 22"`), or from the midpoint of named atoms (`-a "A:315:CA,A:381:CA"`), which is how pockets defined by two residues were set up in the original work.

## Redocking validation, run here

PDB **1IEP** (Abl kinase with imatinib), receptor and crystal ligand from the AutoDock Vina example set, box from the crystal ligand + 6 Å, exhaustiveness 8:

| ligand | best affinity (kcal/mol) | RMSD to the crystal pose |
|---|---|---|
| nilotinib (from SMILES) | −14.17 | different compound |
| imatinib, crystal conformer | −13.22 | **0.40 Å** |
| imatinib, rebuilt from SMILES | −12.04 | **0.90 Å** |
| 2-hydroxybenzophenone (negative control) | −8.60 | different compound |

Redocking the crystal ligand to 0.40 Å, and a structure built from SMILES alone to 0.90 Å, is the standard check that the box, the protonation and the search settings are right: below 2 Å counts as reproducing the experimental binding mode. The unrelated small molecule scores 4.6 kcal/mol worse, as it should. Files in `examples/docking-1iep/`.

## Consensus ranking

When the same ligands are docked with several programs, no single score is reliable, and the ranks agree more than the values do. `consensus-rank.sh` reads one CSV with a column per program and computes four combinations: the average of auto-scaled scores (AASS), the average Z-score, rank by rank (RbR) and exponential consensus ranking (ECR), following Palacio-Rodríguez, Lans, Cavasotto and Cossio, *Sci. Rep.* **2019**, 9, 5142. Programs whose score is better when **high** (GOLD fitness, PLP) are declared with `-H plp:high`; missing scores are skipped for that ligand instead of invalidating the row.

## Corrections made to the original scripts

* **MGLTools** (`prepare_ligand4.py`, Python 2) was replaced by **Meeko**, which is the current preparation route for Vina 1.2 and is maintained.
* **Metal charges** were patched with `sed` on the PDBQT after preparation. That silently depended on the atom order; metal handling now belongs to the receptor preparation step, and the script says so rather than editing files blindly.
* **Missing scores.** The consensus script handled `NIL` entries with a chain of `grep -v NIL` branches that produced different column counts for different ligands. The new one keeps the row and averages over the programs that have a score.
* **Score direction.** The old script assumed every score was better when lower, which silently inverted the ranking of MolDock/PLP style fitness values. The direction is now declared per program.
* **Temporary files.** The original wrote dozens of `ss`, `ss1`, `ss2` files into the working directory and generated small shell scripts on the fly with `echo` and then ran them; everything is now one awk pass.
* **Redocking RMSD** is new: it is the cheapest validation available and it was not in the original workflow.

## Tested

`tests/test-13-docking.sh` reproduces the 1IEP preparation and docking (short exhaustiveness), checks that the crystal ligand redocks within 2 Å, and checks the consensus ranking against a hand-computed example. HADDOCK and Glide/MOE are not installed here, so their helpers are tested on their input and output formats only.
