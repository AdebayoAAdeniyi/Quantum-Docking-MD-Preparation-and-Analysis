# 02 · ORCA 5/6

Ground and excited state workflows for ORCA: optimisation and thermochemistry, TD-DFT with spin-orbit coupling, CASSCF/NEVPT2, fluorescence and phosphorescence rates (ESD module) and NEB-TS reaction paths.

```
02-orca/
├── prepare/
│   ├── orca-input.sh       opt, sp, freq, temps, td, soc, nmr, excited (S0 + T1), casscf (3 steps), neb
│   └── orca-esd-input.sh   ESD(FLUOR) and ESD(PHOSP, 3 sublevels) from finished S0/T1 jobs, DELE computed
├── run/
│   └── orca-status.sh      per step: finished, error, SCF/optimisation convergence, imaginary modes
├── extract/
│   └── orca-extract.sh     CSV: summary, thermo, charges, freq, uv, socuv, socme, casscf, neb, esd, nmr, geometry
└── plot/
    └── plot-neb.sh         energy profile along the NEB path with CI and TS marked
```
IR, Raman, UV-Vis and NMR plots use the shared `bin/qdm-plot-spectra`, which reads the same CSV columns for Gaussian and ORCA.

## Workflow

```bash
# ground state, UV-Vis with and without SOC, NMR, all in methanol (SMD)
02-orca/prepare/orca-input.sh -l list -d xyz -j opt,temps,td,soc,nmr -s smd:methanol
bin/qdm-submit -l list -p orca-jobs/
02-orca/run/orca-status.sh -l list
02-orca/extract/orca-extract.sh -k all -l list
bin/qdm-plot-spectra -t uv -i orca-results/uv.csv

# emission: optimise S0 and T1, then compute the rates
02-orca/prepare/orca-input.sh -l list -d xyz -j excited -s smd:thf
02-orca/prepare/orca-esd-input.sh -l list -s smd:thf        # after the S0/T1 jobs finish
bin/qdm-submit -s run-esd.sh -l list -p orca-jobs/
02-orca/extract/orca-extract.sh -k esd -l list

# multireference: (6,5) active space, 4 singlets + 3 triplets, NEVPT2 with SOC
02-orca/prepare/orca-input.sh -l list -j casscf -a 6,5 -r 4,3 -u 1,3 -S
02-orca/extract/orca-extract.sh -k casscf -l list    # warns about active orbitals with n < 0.02 or > 1.98

# reaction path: NAME-r.xyz and NAME-p.xyz with the same atom order
02-orca/prepare/orca-input.sh -l reactions -d xyz -j neb -i 8
02-orca/extract/orca-extract.sh -k neb -l reactions && 02-orca/plot/plot-neb.sh
```


The complete sequence, a series of molecules through DFT to a statistics table, is in [`../docs/running.md`](../docs/running.md#a-quantum-chemistry-on-a-series-of-molecules). Every script also prints its own usage block with `-h`.

## Design notes

* **One input per step** (`NAME-opt.inp`, `NAME-td.inp` ...). Later steps start from `NAME-opt.xyz` and read `NAME-opt.gbw` with `MORead`, so they are quick, and every output can be parsed on its own without splitting.
* **Memory.** `%maxcore` is set per process to 75 % of the job memory divided by the process count, following the ORCA manual's advice to leave headroom.
* **Heavy atoms.** def2 basis sets include the def2 ECPs from Rb onwards, so complexes need no special basis input (unlike the Gaussian `Gen` route in module 01).
* **Solvation.** `-s smd:SOLVENT` writes `CPCM(SOLVENT)` plus a `%cpcm smd true` block; `-s cpcm:SOLVENT` uses plain CPCM.
* **CASSCF.** Step 1 gives HF orbitals with `NormalPrint` so the candidate active orbitals can be inspected; step 2 converges a state-averaged CASSCF; step 3 re-reads those orbitals and adds SC-NEVPT2 for every requested multiplicity, with optional QDPT spin-orbit coupling. The extractor prints the natural occupations and flags orbitals that are barely correlated, which was a manual `grep 'N(occ)'` check in the original notes.
* **Phosphorescence.** ORCA needs the adiabatic S0–T1 energy (`DELE`) as input and computes each triplet sublevel separately ([ORCA manual, Excited State Dynamics](https://www.faccts.de/docs/orca/6.0/manual/contents/typical/esd.html)). `orca-esd-input.sh` reads both energies, converts the difference to cm⁻¹ and writes the three `IROOT` jobs.

## Corrections made to the original scripts

* The NEB, CASSCF, SOC, fluorescence and phosphorescence generators were separate scripts with the molecule names typed inside them; they are now steps of one generator that reads a list.
* `%maxcore 10000` with 24 processes asked for 240 GB per node regardless of the machine; memory now comes from the cluster settings.
* The S0–T1 energy for `DELE` was taken from "Electronic energy" plus "Total thermal correction", which adds the 298 K thermal energy of each state to an adiabatic gap. It now uses the electronic energies, and `-Z` adds the zero point energies (the usual adiabatic 0–0 definition) when wanted.
* TD-DFT tables changed layout in ORCA 6 (states are labelled `0-1Ag -> 1-1Bu`, energies in eV first). The parser reads both layouts and matches every state to its orbital transitions by energy.
* `neb_snapshots.py` (third party) is not included.

## Tested

`tests/test-02-orca.sh` (27 checks): inputs for the five benzophenones; parsers on public ORCA 4.2, 5.0 and 6.0 outputs from [cclib](https://github.com/cclib/cclib) (thermochemistry, both TD-DFT layouts, IR and Raman, NMR, Hirshfeld charges of a radical, spin-orbit matrix elements and SOC corrected spectrum, CASSCF/NEVPT2 transition energies, an SCF failure). No public NEB or ESD output was available, so those two parsers are tested on short fragments that copy the layout printed in the ORCA tutorial and manual (`tests/fixtures/`). ORCA itself was not run.
