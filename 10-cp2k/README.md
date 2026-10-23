# 10 · CP2K: classical and QM/MM molecular dynamics

CP2K runs MD with the same Amber topology built in module 08, either fully classical (FIST) or with part of the system treated by DFT and the rest by the force field (QM/MM with GEEP electrostatic embedding). That makes it the natural next step when a reaction, a metal centre or a proton transfer has to be described inside a solvated system that was equilibrated with Amber.

```
10-cp2k/
└── cp2k-input.sh    -j mm    classical MD (FIST, SPME, NPT_I or NVT)
                     -j qmmm  QM/MM MD: PBE-D3/DZVP-MOLOPT-SR (GPW) inside the force field
```

## Workflow

```bash
# after module 08 has built and equilibrated the system
10-cp2k/cp2k-input.sh -d amber-systems/complex -j mm -n 50000 -E NPT_I
bin/qdm-submit -s run-cp2k-mm.sh amber-systems/complex

# QM/MM: the ligand (residue LIG) is the QM region, charge +1
10-cp2k/cp2k-input.sh -d amber-systems/complex -j qmmm -q LIG -Q 1 -n 20000 -t 0.5
# metadynamics on top (module 12)
10-cp2k/cp2k-input.sh -d amber-systems/complex -j qmmm -q LIG -P
```

The QM region is given either as a residue name (`-q LIG`) or as atom indices (`-q "1-25,40"`). The script finds the atoms in the prmtop itself.

## Things this script handles that cost real debugging time

* **"Unknown element for KIND"**. CP2K builds its atomic kinds from the Amber atom names and force field types, and neither is an element symbol: GAFF2 writes types such as `nu`, `hx` and `ca`, and atom names such as `C12`. The script writes an explicit `&KIND` with `ELEMENT` for every distinct name **and** type in the topology. This was the first error hit when running the original CP2K input on a current system.
* **NetCDF restarts.** Amber 20+ writes binary restarts, which CP2K cannot read. The coordinates are converted with cpptraj to a PDB, which also carries the element column and the box as `CRYST1`.
* **QM cell.** It is sized from the actual extent of the QM region plus 6 Å on each side, and kept cubic because the wavelet Poisson solver for non-periodic QM cells requires it.
* **Basis sets for QM kinds.** The kinds that contain QM atoms, and the kind named after each QM element, get `DZVP-MOLOPT-SR-GTH` and `GTH-PBE`.
* **Cut bonds.** If the QM selection does not cover whole residues, the script warns that `&LINK` sections are needed; CP2K otherwise stops with a connectivity error.
* `EI_SCALE14 0.8333333` and `VDW_SCALE14 0.5` are set, which are the Amber 1-4 scalings. Leaving them at the CP2K defaults silently changes the force field.

## Tested

Both paths were run here on a real solvated system (imatinib in water, 5,363 atoms, built in module 08):

* **MM**: 20 steps of NVT with the Amber force field, SPME electrostatics; total FIST energy −28.5 Eh, temperature around 290 K.
* **QM/MM**: 2 steps with one water molecule as the QM region (PBE-D3, DZVP-MOLOPT-SR, 400 Ry, GEEP embedding, wavelet Poisson solver); total QM/MM energy −45.4 Eh.

A 69-atom QM region was also set up correctly but ran out of memory in this 7 GB container, which is expected: a QM region that size needs a proper node. `tests/test-10-cp2k.sh` runs the MM case and checks the QM/MM input.

The ASE slab-building scripts that sat beside the original CP2K inputs used the old `ase.lattice.surface` and `ase.structure` modules, which were removed from ASE years ago. They are not carried over; `ase.build` is the current route and is a few lines in a notebook rather than a workflow script.
