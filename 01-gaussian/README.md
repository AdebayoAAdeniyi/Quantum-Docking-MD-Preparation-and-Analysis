# 01 · Gaussian 09/16

Batch preparation, submission, extraction and plotting for Gaussian jobs. A whole series of molecules is handled in one command: every molecule gets its own folder, one input with all the requested steps, and a job script for a laptop, a PBS cluster or a SLURM cluster.

```
01-gaussian/
├── prepare/
│   ├── gaussian-input.sh       multi-step jobs (opt+freq, NPA/Hirshfeld, CHelpG, TD-DFT, NMR, polar, EPR, wfx) + Fukui inputs
│   ├── gaussian-ts-input.sh    TS search + IRC, QST3, counterpoise (BSSE) jobs
│   └── gaussian-followup.sh    restarts from a checkpoint: thermochemistry at many temperatures,
│                               TD-DFT in several solvents, NBO/NRT, wavefunction files
├── run/
│   └── gaussian-status.sh      finished / running / failed, imaginary frequencies, restart geometries
└── extract/
    ├── split-link1-log.sh      cut a --Link1-- log into one file per step
    └── gaussian-extract.sh     CSV tables: summary, thermo, charges, uv, freq, nmr, jcoupling,
                                polar, hyperfine, fukui, geometry
```
Plots are made with the shared [`bin/qdm-plot-spectra`](../bin/qdm-plot-spectra) (IR, Raman, UV-Vis, NMR).

## Workflow

```bash
# 1. structures in any format Open Babel reads, one name per line in "list"
cat list
#   1-H
#   2-OMe
#   MoOCl4
cat chargemult.txt            # only molecules that are not neutral singlets
#   MoOCl4 -1 2

# 2. inputs: optimisation + frequencies in the gas phase, properties in acetonitrile,
#    plus the N, N+1, N-1 single points for conceptual DFT
01-gaussian/prepare/gaussian-input.sh -l list -d structures \
    -j opt,pop,chelpg,td,nmr -m B3LYP -b "6-311+G(d,p)" -e LanL2DZ -s smd:acetonitrile -F

# 3. run (scheduler, cores, memory and modules come from config/cluster.conf)
bin/qdm-submit -l list -p gaussian-jobs/
01-gaussian/run/gaussian-status.sh -l list -r     # check; failed jobs -> restart/NAME.xyz

# 4. extract everything to CSV and plot
01-gaussian/extract/gaussian-extract.sh -k all -l list -R tms/tms.log
bin/qdm-plot-spectra -t ir -i gaussian-results/freq.csv -s 0.967 -e exp/1-H-ir.dat
bin/qdm-plot-spectra -t uv -i gaussian-results/uv.csv
bin/qdm-plot-spectra -t nmr -i gaussian-results/nmr.csv -E C
```

## What each step means

| step | route keywords | gives |
|---|---|---|
| `opt` / `optraman` | `Opt Freq` (`Freq=Raman`) | minimum, IR (and Raman), ZPE, H, G, S, Cv |
| `sp` | single point | energy at the input geometry |
| `pop` | `Pop=(Full,NPA,Hirshfeld)` | NPA, Hirshfeld and CM5 charges, all orbitals |
| `chelpg` | `Pop=CHelpG` | charges fitted to the electrostatic potential (MEP) |
| `td` | `TD=(NStates=20)` | vertical excitations, oscillator strengths, orbital contributions |
| `nmr` | `NMR=GIAO` (`-J` adds spin-spin) | shielding tensors, J couplings with FC/SD/PSO/DSO terms |
| `polar` | `Polar` | static α and β (nonlinear optics) |
| `epr` | `Prop=EPR` | hyperfine couplings for radicals |
| `wfx` | `Output=WFX` | wavefunction for Multiwfn and AIMAll (module 05) |
| `-F` | three `Pop` single points | condensed Fukui functions, IP, EA, μ, η, ω |

All steps after the first read geometry and orbitals from the shared checkpoint. When heavy atoms are present (transition metals, lanthanides, Sn, Sb, I, Pb ...), the first step writes a `Gen` basis with the ECP, and later steps use `ChkBasis`, which according to the Gaussian manual also retrieves the ECP from the checkpoint. That keeps one copy of the basis set per input.

## Conceptual DFT definitions used by `-k fukui`

With E(N), E(N+1) and E(N−1) at the same geometry:

* IP = E(N−1) − E(N), EA = E(N) − E(N+1)
* chemical potential μ = −(IP + EA)/2, hardness η = (IP − EA)/2, softness S = 1/(2η) (Parr and Pearson, *J. Am. Chem. Soc.* 1983, 105, 7512)
* electrophilicity ω = μ²/(2η) (Parr, von Szentpály and Liu, *J. Am. Chem. Soc.* 1999, 121, 1922)
* condensed Fukui functions from atomic charges q: f⁺ = q(N) − q(N+1), f⁻ = q(N−1) − q(N), f⁰ = (q(N−1) − q(N+1))/2 (Yang and Mortier, *J. Am. Chem. Soc.* 1986, 108, 5708), and the dual descriptor f⁺ − f⁻ (Morell, Grand and Toro-Labbé, *J. Phys. Chem. A* 2005, 109, 205). Both NPA and Hirshfeld values are reported; Hirshfeld based indices are usually preferred because they are rarely negative.

## Corrections made to the original scripts

* **Global reactivity descriptors.** The old Fukui extraction subtracted the energies in the wrong order, so IP and EA came out with the wrong sign, the column titles called the N+1 species a cation, and the electrophilicity was computed as μ²/2 × η instead of μ²/(2η). All three are fixed above.
* **NMR element selection.** `grep C` also matched Cl and `grep H` matched every line with an H in it; the parser now uses the element column.
* **Mixed steps.** "Last match in the log" mixed gas phase and solvent values from different Link1 steps. Every log is now split into steps and each row says which step it came from.
* **Open shell NPA.** For radicals NBO prints total, α and β summaries; the last one (β only) used to be taken. Only the total density summary is used now.
* **Thermochemistry at several temperatures** now uses `Freq=ReadFC`, which re-uses the stored force constants instead of recomputing the Hessian for each temperature.
* **NMR references.** The old reference shieldings (C 191.12, H 31.42, N −128.97 ppm and the scaled forms) are only valid for the level of theory they were computed at. They can still be given with `-r`, but the recommended way is `-R tms.log`, a reference molecule computed at the same level.
* **β (first hyperpolarisability).** β_tot = (β_x² + β_y² + β_z²)^½ with β_i = β_iii + β_ijj + β_ikk, from the input orientation block, 1 au = 8.6393 × 10⁻³³ esu. Gaussian's own `||` value is printed next to it, since several conventions exist in the literature and papers should say which one they use.
* Hardcoded home and cluster paths, a separate copy of the transpose program in every script, and per-element `grep` atom reordering were replaced by the shared library.

## Tested

`tests/test-01-gaussian.sh` (29 checks) builds inputs for the five benzophenones in `examples/` and a Mo oxo complex, and checks every parser against public Gaussian 16/09/03 outputs from the [cclib](https://github.com/cclib/cclib) test set (thermochemistry, TD-DFT, NMR and J couplings, α and β, Fermi contact couplings, Hirshfeld/CM5/NPA charges, CCSD(T) energies, an unconverged optimisation). Gaussian itself is commercial and was not run; the values in the tests are read directly from those reference outputs.
