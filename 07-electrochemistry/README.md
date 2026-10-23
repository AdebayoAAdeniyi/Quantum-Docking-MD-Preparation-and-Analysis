# 07 · Electrochemistry: computed and measured

Two halves of the same question. `computational/` turns quantum chemistry output into reduction potentials and pKa values through a thermodynamic (Born-Haber) cycle; `experimental/` processes the cyclic voltammetry and UV-Vis files from the bench and gives frontier orbital energies that can be compared with the computed ones.

```
07-electrochemistry/
├── computational/
│   ├── redox-table.sh   builds energies.csv (G_gas, E_gas, E_solv) from any summary.csv
│   └── redox-pka.sh     reduction potentials vs SHE/SCE/Fc, and pKa values
└── experimental/
    ├── cv-process.sh    ferrocene referencing, onset by tangents, E_HOMO / E_LUMO
    ├── uv-onset.sh      lambda_max, absorption onset, optical gap
    └── LAB-PROTOCOL.md  how the CV and UV data were measured
```

## Computed reduction potentials and pKa

The cycle needs three numbers per species: the gas phase free energy (optimisation + frequencies), and single point energies in the gas phase and in implicit solvent at the same geometry.

```bash
# 1. gas phase opt+freq, then single points in gas and in solvent (modules 01, 02 or 04)
01-gaussian/prepare/gaussian-input.sh -l species -j opt                      -o gas
01-gaussian/prepare/gaussian-input.sh -l species -j sp  -m M06 -b "6-311+G(2df,2p)" -o gasSP
01-gaussian/prepare/gaussian-input.sh -l species -j sp  -m M06 -b "6-311+G(2df,2p)" -s smd:acetonitrile -o solv
# 2. one energy table
07-electrochemistry/computational/redox-table.sh -g gas/summary.csv -H gasSP/summary.csv -s solv/summary.csv > energies.csv
# 3. potentials (O + n e- -> R) and pKa (HA -> A- + H+)
cat reactions.txt
#   quinone  BQ  BQ-  1
07-electrochemistry/computational/redox-pka.sh -k redox -t energies.csv -r reactions.txt -R Fc
07-electrochemistry/computational/redox-pka.sh -k pka   -t energies.csv -r acids.txt
```

Constants, all citable and all overridable from the environment:

| quantity | value | source |
|---|---|---|
| absolute potential of the SHE | 4.28 V | Isse and Gennaro, *J. Phys. Chem. B* 2010, 114, 7894 |
| E(SCE) vs SHE | +0.241 V | standard tables |
| E(Fc⁺/Fc) vs SHE in acetonitrile | +0.63 V | Pavlishchuk and Addison, *Inorg. Chim. Acta* 2000, 298, 97 |
| G_gas(H⁺) | −6.28 kcal/mol | Sackur-Tetrode, 298 K |
| ΔG_solv(H⁺) in water | −265.9 kcal/mol | Tissandier *et al.* 1998; Kelly, Cramer and Truhlar 2006 |
| 1 atm → 1 M | RT ln(RT/p°) = 1.89 kcal/mol at 298 K | standard state conversion |

## Experimental CV and UV

```bash
07-electrochemistry/experimental/cv-process.sh -F FC.TXT -d red -g 3.6 -w 200:700 1-H.TXT 2-OMe.TXT
07-electrochemistry/experimental/uv-onset.sh -w 250:600 *.txt
```

* `cv-process.sh` reads the potentiostat text export whatever the header looks like (numeric lines only), splits the scan into cycles and takes the last complete one, locates the ferrocene anodic and cathodic peaks in the window given with `-w`, reports E½ and ΔEp (with a warning when ΔEp is far from the 59 mV of a reversible one-electron couple), shifts every analyte to the Fc⁺/Fc scale, and finds the onset as the intersection of the baseline tangent with the tangent at the steepest point of the rising current.
* E_LUMO = −(E_red,onset + 4.8) eV and E_HOMO = −(E_ox,onset + 4.8) eV; with `-g Eg` the other orbital follows from the optical gap. The 4.8 eV value for Fc⁺/Fc below vacuum is the common choice but not the only one, so `-E` changes it and the value used is written in the output.

## Corrections made to the original scripts

* **Hardcoded windows.** The old scripts contained tables of peak-search windows and tangent ranges for specific compounds, so they could not be reused. Peaks and the steepest point are now found automatically inside a window given on the command line.
* **Cycle selection.** The old route cut the file with `csplit` on a value copied from the data by hand. Cycles are now detected from the sweep direction, and an incomplete last cycle falls back to the previous one.
* **Fixed header length.** `tail -n +26` broke whenever the instrument wrote a different header. Any non-numeric line is now skipped.
* **A check that was missing.** ΔEp of the ferrocene couple is reported, which tells you at once whether the reference behaved reversibly in that solvent and electrolyte.
* **pKa sign conventions.** The old notes said to multiply by −1 in the deprotonation case; the reaction is now always written HA → A⁻ + H⁺, so no sign juggling is needed.

## Tested

`tests/test-07-electrochemistry.sh`: the CV and UV processing is checked against simulated data whose true values are known by construction (a reversible couple with E½ exactly 0.450 V and peaks ±35 mV, and a reduction wave whose tangent onset is analytically E₀ + 2s); the Born-Haber arithmetic is checked against a hand calculation. **No real potentiostat files were available**, so the experimental scripts have not been run on instrument data; the file reader is written to tolerate header blocks, `Segment` lines, CRLF endings and comma decimal separators, but it should be checked on your own export once.
