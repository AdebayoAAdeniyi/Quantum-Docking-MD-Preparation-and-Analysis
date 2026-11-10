# 04 · NWChem and xTB (free programs, run for real)

**NWChem** handles the open-source DFT work in the collection: optimisations with thermochemistry over a range of temperatures, transition state searches, NEB reaction paths and counterpoise-corrected interaction energies. **GFN2-xTB** is the fast screening layer. It pre-optimises structures before DFT and docking, and it ranks a series by gap, IP, EA and reactivity indices in seconds per molecule. That replaces the old practice of running quick PM3 pre-optimisations through Gaussian under wine.

```
04-nwchem-xtb/
├── prepare/nwchem-input.sh     opt (+freq at many T), sp, ts (freq -> saddle -> freq), neb, bsse
├── run/
│   ├── xtb-run.sh              opt / ohess / vipea / vfukui for a whole set, in seconds
│   └── nwchem-displace.sh      geometries displaced along an imaginary mode (restart or TS check)
└── extract/
    ├── nwchem-extract.sh       CSV: summary, thermo, freq, neb, bsse, geometry
    └── xtb-extract.sh          one screening table (+ mu, eta, omega) and a Fukui table
```

## Workflow

```bash
# 1. screen everything with xTB in implicit water
04-nwchem-xtb/run/xtb-run.sh -s water -l list -d structures
04-nwchem-xtb/extract/xtb-extract.sh            # xtb-results/xtb-screen.csv, xtb-fukui.csv

# 2. DFT on the selected molecules, starting from the xTB geometries
04-nwchem-xtb/prepare/nwchem-input.sh -j opt -b "6-311+g(2d,p)" -D -d xtb-jobs-geoms -l selected
bin/qdm-submit -s run-opt.sh -l selected -p nwchem-jobs/

# 3. a reaction: TS search from NAME-ts, check it, then kinetics (module 06)
04-nwchem-xtb/prepare/nwchem-input.sh -j ts -l reactions -d xyz
04-nwchem-xtb/run/nwchem-displace.sh nwchem-jobs/rx/rx-ts.log      # rx-plus.xyz / rx-minus.xyz
04-nwchem-xtb/extract/nwchem-extract.sh -k all
```


The complete sequence, screening a series with xTB and taking the survivors to DFT, is in [`../docs/running.md`](../docs/running.md#a-quantum-chemistry-on-a-series-of-molecules). Every script also prints its own usage block with `-h`.

## Results produced while building this module

**xTB, five 2-hydroxybenzophenones** (UV absorber scaffolds built from SMILES in `examples/`, GFN2-xTB with ALPB water, `examples/xtb-benzophenones/`):

| molecule | gap (eV) | dipole (D) | IP (eV) | EA (eV) | ω (eV) |
|---|---|---|---|---|---|
| 1-H | 2.869 | 6.56 | 7.667 | 3.001 | 6.10 |
| 2-OMe | 2.722 | 6.97 | 7.569 | 3.056 | 6.25 |
| 3-Oct | 2.719 | 7.00 | 7.490 | 3.045 | 6.24 |
| 4-Cl | 2.747 | 7.65 | 7.448 | 3.126 | 6.47 |
| 5-Br | 2.758 | 8.88 | 7.516 | 3.056 | 6.27 |

All five are minima (no imaginary modes). The xTB "gap" is a semi-empirical orbital gap and is used only for ranking, not as an optical gap.

**NWChem, B3LYP/6-31G\*** (`examples/kinetics-hcn-hnc/`):
* Water optimised with frequencies at 298.15 and 310 K: E = −76.40702 Eh, three real modes.
* The HCN → HNC transition state was found with `task dft saddle` from a bent guess. It has one imaginary mode at −1147 cm⁻¹, and the ±100 % displaced structures point to the right minima: one has H–C 1.17 Å (towards HCN), the other H–N 1.06 Å (towards HNC).
* The water dimer counterpoise job gives a BSSE of 0.87 kcal/mol at this level.
* The HNC − HCN energy is 16.1 kcal/mol including ZPE. The high-level reference is 5212 ± 30 cm⁻¹ = 14.9 kcal/mol at 0 K (*J. Phys. Chem. A* 2015, 119, 10929), so B3LYP/6-31G\* is about 1.2 kcal/mol too high. That is the expected size of error for this level and a useful reminder to benchmark.
* **Not converged:** an NWChem NEB on the same reaction (8 beads, starting path through the TS guess) ran 50 iterations. The bead next to the TS reached 45 kcal/mol, but the beads near the end points were still far too high, so that path was not converged. The saddle search is the reliable route for this kind of reaction, and the README says so rather than showing an unconverged path.

## Corrections made to the original scripts

* **Molecule names inside scripts.** Each old NWChem script had the molecule names, charges and multiplicities typed into it. They now come from a list and `chargemult.txt`.
* **NEB starting path.** With a straight-line interpolation, the moving H atom passed through the C and N atoms, and the first NEB run started at about 1000 kcal/mol. When a `NAME-ts` file is present, the starting path is now built in two halves, R → TS → P, and passed with `xyz_path`. NWChem reads that file from `permanent_dir`, which the script takes care of.
* **Counterpoise ghost atoms.** The counterpoise input failed with "no basis on ghost atom". Ghost atoms (`bqX`) now get the same basis as the real atoms.
* **TS search Hessian.** `task dft saddle` now reads the Hessian from the preceding frequency job (`inhess 2`), so the search follows the right mode from the start.
* **Gibbs energy.** NWChem prints H and S but not G. The extractor computes G = E + H_corr − TS, and reports the 298.15 K values when that temperature was computed.
* **Cluster details.** Cluster paths, project codes and modules are in `config/cluster.conf`, not in each script.

## Tested

`tests/test-04-nwchem-xtb.sh` runs xTB on two benzophenones and an NWChem water optimisation with frequencies (10 checks). `QDM_SLOW=1` adds the HCN/HNC transition state search and the displacement check.
