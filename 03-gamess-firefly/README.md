# 03 · GAMESS-US and Firefly: NBO energy decomposition (NEDA)

I used GAMESS and Firefly mainly for one thing that Gaussian cannot do: the **natural energy decomposition analysis (NEDA)** of NBO. It splits the interaction between a metal and its ligands, or between two molecules, into electrical (electrostatic + polarisation), charge transfer and core (exchange + deformation) terms. NEDA runs only inside GAMESS, Firefly and NWChem ([NBO tutorial](https://nbo6.chem.wisc.edu/tut_neda.htm)).

```
03-gamess-firefly/
├── prepare/gamess-input.sh    energy, opt, hessian, optfreq, neda; GAMESS-US (-P gamess) or Firefly (-P firefly)
└── extract/gamess-extract.sh  CSV: summary, thermo, freq, charges, neda (+ neda-wide.csv), geometry
```

## Workflow

```bash
pip install basis_set_exchange          # basis sets and ECPs are taken from the BSE library

# fragments given as atom ranges ...
03-gamess-firefly/prepare/gamess-input.sh -P firefly -j neda -f "1-12/13-30" -l list -d xyz
# ... or as separate files NAME-f1.pdb, NAME-f2.pdb cut from the complex in any viewer;
# their atoms are matched to the complex by coordinates
03-gamess-firefly/prepare/gamess-input.sh -j neda -l list -d structures -b "6-31G*" -e LANL2DZ

bin/qdm-submit -s run-neda.sh -l list -p gamess-jobs/
03-gamess-firefly/extract/gamess-extract.sh -k neda -l list     # neda.csv and neda-wide.csv
```

## How the inputs are built

* Every atom carries its own basis functions in `$DATA` and every heavy atom its potential in `$ECP` (`PP=READ` for GAMESS-US, `ECP=READ` for Firefly), with `NONE` for light atoms and the name alone for a repeated potential. Basis and ECP data come from the [Basis Set Exchange](https://www.basissetexchange.org), so any combination such as 6-31G\* on C, H, N, O, Cl and LANL2DZ on Ru, Pt or I can be used. The original workflow kept a folder of hand-copied library files and spliced them in with `sed`.
* For an open-shell system the NEDA request is repeated for the α and β spins, as the NBO manual requires.
* `$SYSTEM MWORDS` is set from the job memory (words of 8 bytes per process, with 25 % headroom).

## Corrections made to the original scripts

* **NEDA fragments.** Some of the Firefly inputs wrote `NEDA (1,2,...) (...)`. Without the `NAO` keyword NBO reads those numbers as *molecular units*, not atoms. The atom lists now always go with `NAO NEDA`, which is the documented form for atom-based fragments.
* The NBO licence key that was stored in the old Firefly header files is not part of this repository. NBO must be licensed and linked separately.
* The old NEDA extraction searched for the label `XC`, which older NBO versions printed; NBO 6 prints `EX`. Both are read.

## Tested

`tests/test-03-gamess.sh` (18 checks). The Mo basis set and ECP generated for MoOCl₄ are identical to the ones in the working GAMESS-US input from the [cclib](https://github.com/cclib/cclib) test set. The parsers are checked on cclib GAMESS-US 2017 and Firefly 8 outputs (energies, ZPE, HOMO, a converged optimisation, an SCF failure, an NH₃ transition state with one imaginary mode, a dianion), and the NEDA parser on the water dimer table printed in the NBO 6 tutorial. GAMESS, Firefly and NBO were not run here.
