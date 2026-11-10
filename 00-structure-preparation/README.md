# 00 · Structure preparation

Everything downstream starts here. A quantum chemistry job needs clean 3D
coordinates with the right protonation, `antechamber` and `tleap` need a ligand
PDB with unique atom names, docking needs one file per molecule, and all of
them need whatever format you were given turned into whatever format the next
program reads.

These four scripts do that, through Open Babel.

```
00-structure-preparation/
├── convert-structures.sh     any format to any other, hydrogens, pH, centring
├── build-3d-from-smiles.sh   a SMILES table to optimised 3D conformers
├── ligand-pdb-for-md.sh      unique atom names and one residue name, for MD
└── split-multistructure.sh   a multi-molecule file to one file per molecule
```

## Workflow

```bash
# a series of molecules from SMILES: protonate for pH 7.4, 30 conformers,
# MMFF94 minimisation, written as xyz, pdb and mol2
00-structure-preparation/build-3d-from-smiles.sh -i smiles.txt -o prepared \
    -f xyz,pdb,mol2 -p 7.4 -n 30

# convert a folder of crystal structures, adding hydrogens at pH 7.4
00-structure-preparation/convert-structures.sh -l list -d raw -o prepared \
    -f xyz,pdb -p 7.4

# a docked pose into a ligand PDB that antechamber and tleap will accept
00-structure-preparation/ligand-pdb-for-md.sh best-pose.sdf LIG > ligand.pdb

# split a ZINC download, a docking output or an NMR ensemble
00-structure-preparation/split-multistructure.sh poses.sdf -o split -k 1
```

The SMILES table is two columns, name and SMILES:

```text
1-H    Oc1ccccc1C(=O)c1ccccc1
2-OMe  COc1ccc(C(=O)c2ccccc2)c(O)c1
```

The name becomes the file name and the key in `list` and `chargemult.txt`, so
it carries through every later module. Output is flat, one file per molecule
per format: `built/benzene.xyz`, `built/benzene.pdb` and so on, with the
force-field energies of the retained conformers in `built/build-energies.txt`.

## What each script does, and why it matters later

**`convert-structures.sh`** reads everything Open Babel reads, including
Gaussian `.log` and `.com` files, so the optimised geometry of one calculation
becomes the input of the next without a manual step. `-p pH` protonates for a
pH rather than just adding hydrogens, which is the difference between a
carboxylic acid and a carboxylate and therefore between a sensible docking
result and a meaningless one. `-c` centres the molecule at the origin.

**`build-3d-from-smiles.sh`** generates 3D coordinates, searches conformers with
a genetic algorithm, and minimises with MMFF94, falling back to UFF for
elements MMFF94 does not cover. It writes the lowest-energy conformer in every
requested format. `-n` sets the number of conformers and `-s` the minimisation
steps.

A SMILES string carries no 3D information, so the conformer this produces is a
reasonable starting point and nothing more. For anything where the conformer
matters, re-optimise it with xTB (module 04) or DFT (modules 01, 02, 04) before
drawing conclusions.

**`ligand-pdb-for-md.sh`** gives every atom a unique name of the form C1, C2,
O1, H1 and puts the whole molecule in one residue. This is not cosmetic:
`antechamber`, `tleap` and `pdb2gmx` identify atoms by name, and a PDB with two
atoms called `C` or with the ligand split across several residue numbers will
either fail or, worse, produce a topology that is subtly wrong. It converts the
input first, so an sdf or mol2 works as well as a pdb.

**`split-multistructure.sh`** names the output files after the molecule title,
such as a ZINC identifier, or the MODEL number, and appends `_01`, `_02` to
repeated titles so no pose is overwritten. `-k 1` keeps only the first
structure per title, which is the best pose in a docking output.

## Corrections made to the original scripts

* The original route to 3D structures, `download-2d-nist-optimized-to-3d-obabel`,
  fetched 2D structures from the NIST WebBook and then ran three `obminimize`
  passes by hand. It needed internet access to run at all.
  `build-3d-from-smiles.sh` works from a SMILES table offline and does the
  conformer search in one call.
* Four separate splitters, `splitsdf.csh`, `splitmol2.csh`, `splitpdb.csh` and
  the MOE and Glide splitters, each did the same job for one format with
  `csplit`. They are one script that handles every format.
* `numbering-pdb` and `numbering-2-pdb` renamed ligand atoms for MD but
  overwrote repeated names, which silently merged atoms in the topology.
  `ligand-pdb-for-md.sh` guarantees uniqueness per element.
* Protonation was done by adding hydrogens without reference to pH, which is
  correct only for neutral molecules. The pH is now an argument.

## Tested

`tests/test-00-structure-preparation.sh`, 18 checks, all running real Open
Babel calls:

* two molecules built from SMILES, with the atom counts verified against the
  molecular formulas and the coordinates confirmed to be 3D rather than flat;
* a format round trip that checks the atom count and the element counts survive;
* centring, checked by computing the centroid;
* pH protonation, checked by confirming acetic acid carries fewer hydrogens at
  pH 11 than at pH 2, which is the behaviour that matters for docking;
* the MD ligand naming, checked for atom count, name uniqueness, a single
  residue name and number, and names that begin with their element symbol;
* splitting, checked for nothing lost, a repeated title not overwriting itself,
  and `-k 1` keeping one structure per title.

## Next

Module 04 for fast xTB screening of the series, modules 01, 02 or 04 for DFT,
module 13 for docking, or module 08 to build an MD system. The complete
sequences are in [`../docs/running.md`](../docs/running.md).
