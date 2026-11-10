# 05 · Wavefunction analysis: Multiwfn and QTAIM (AIMAll)

Analyses that start from a finished wavefunction (`.fchk`, `.wfx`) rather than from a new quantum chemistry run: Fukui function densities, the electrostatic potential on the molecular surface, frontier orbital cubes and overlap, hole–electron analysis of excited states, and bond critical point properties from the quantum theory of atoms in molecules (QTAIM).

```
05-wavefunction-analysis/
├── run/
│   ├── multiwfn-batch.sh   fukui | mesp | homo-lumo | overlap | hole-electron | custom key file
│   └── aimall-run.sh       AIMAll without the GUI for every wavefunction in a batch
└── extract/
    └── qtaim-extract.sh    bond critical points: rho, Laplacian, ellipticity, G, V, H, -G/V, bond class
```

## Workflow

```bash
# wavefunctions from module 01: gaussian-input.sh -F (Fukui fchk files), -j ...,wfx
05-wavefunction-analysis/run/multiwfn-batch.sh -t fukui -l list        # f+ and f- density cubes
05-wavefunction-analysis/run/multiwfn-batch.sh -t mesp -l list         # mesp.csv: Vmin, Vmax per molecule
05-wavefunction-analysis/run/multiwfn-batch.sh -t overlap -l list      # HOMO/LUMO centroid distance, overlap
05-wavefunction-analysis/run/multiwfn-batch.sh -t hole-electron -s 1 -l list

05-wavefunction-analysis/run/aimall-run.sh -l list -s -wfx -q         # job scripts; submit with qdm-submit -s run-aimall.sh
05-wavefunction-analysis/extract/qtaim-extract.sh -l list -e Fe,Ru     # bcp-all.csv, bcp-weak.csv, bcp-metal.csv
```

## Interpreting the QTAIM table

At each bond critical point the extractor adds the total energy density H = G + V and the ratio −G/V, and labels the interaction from the sign of the Laplacian ∇²ρ and of H (Cremer and Kraka, *Angew. Chem. Int. Ed.* 1984, 23, 627; Espinosa, Alkorta, Elguero and Molins, *J. Chem. Phys.* 2002, 117, 5529):

| class | ∇²ρ | H | typical |
|---|---|---|---|
| shared | < 0 | < 0 | covalent bonds |
| intermediate | > 0 | < 0 | metal–ligand, dative and strong hydrogen bonds |
| closed-shell | > 0 | > 0 | ionic contacts, ordinary hydrogen bonds, van der Waals |

`bcp-weak.csv` keeps the critical points with ρ ≤ 0.04 a.u. (hydrogen-bond-like contacts) and `bcp-metal.csv` those that involve the elements given with `-e`; these are the two selections the original `select-bonds-qtaim` made.

## Changes from the original scripts

* About twenty `run-Multifwx-*` scripts, many of them copies for three different clusters, are one driver with a task switch. The menu key sequences are the ones those scripts used (Multiwfn 3.8); `-t custom -k keys.txt` runs any other sequence, with `NAME` replaced by each molecule.
* The HOMO index is read from Multiwfn's own start-up screen for every molecule, as before, but the two message formats (closed and open shell) are handled in one place.
* The Fukui cubes now read the `NAME-fukui-N/Nplus1/Nminus1.fchk` files that `gaussian-input.sh -F` produces, so the chain from input to density difference runs without renaming files.
* The BCP extraction used `split -l 8` on grep output, which breaks when AIMAll prints a different number of lines per critical point; the new parser reads each critical point block by its keys.
* Tian Lu's VMD scripts (`orbitalshow.vmd`) and the large saved VMD states are not redistributed; they ship with the Multiwfn manual and example package from [sobereva.com/multiwfn](http://sobereva.com/multiwfn). See [`../docs/third-party-material.md`](../docs/third-party-material.md#tian-lus-vmd-scripts) for how the cube and the rendering step fit together.

## Tested

Multiwfn and AIMAll were not available in the build environment, so this is the least tested module. `tests/test-05-wavefunction.sh` checks the QTAIM classification and filters on a small fragment that copies the key = value layout the original script searched for (values invented), and checks the syntax of every script. The Multiwfn key sequences are taken unchanged from scripts that produced published results, but they should be checked once against your Multiwfn version.
