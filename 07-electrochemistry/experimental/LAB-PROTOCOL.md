# Cyclic voltammetry: bench protocol behind the data

These are the notes I kept for measuring organic and organometallic compounds by cyclic voltammetry (CV), square wave voltammetry (SWV) and linear sweep voltammetry (LSV). They describe how the files that `cv-process.sh` reads were produced, and they have been tidied for other users.

## Cell and solution

1. **Electrolyte.** Use about 200 equivalents of supporting electrolyte relative to the analyte. With 0.5 mM analyte that is 0.1 M, for example tetrabutylammonium hexafluorophosphate (Bu₄NPF₆). A weakly coordinating borate such as [B(C₆F₅)₄]⁻ is an option for very reactive cations. *(The original note said "like 1 M"; 0.1 M is what the 200:1 ratio gives and is the usual concentration.)*
2. **Solvent.** Use about 1 to 10 mL of dry electrochemistry-grade solvent, drawn through a septum with a nitrogen-flushed needle so that no air enters. Dichloromethane or acetonitrile behave better than DMSO.
3. **Electrodes.** Use glassy carbon as the working electrode (a larger surface gives clearer peaks) and platinum wire as the counter and pseudo-reference electrodes.
4. **Polishing.** Rinse the glassy carbon with methanol and acetone and dry only its sides. Then polish the tip on a wet pad with alumina paste, rinse with methanol and acetone again, and do not touch the surface.
5. **Degassing.** Bubble nitrogen or argon through the solution before every series and keep a blanket over it during the scans. Dissolved O₂ gives its own reduction wave.

## Measurement sequence

1. **Blank (solvent + electrolyte).** Find the potential window and set the scan limits inside it.
2. **Analyte (0.5 mM).** Degas, then record CVs at 100, 200, 300, 400, 500 and 1000 mV s⁻¹.
3. **Internal reference.** Add ferrocene at the same concentration as the analyte (0.5 mM; 0.093 mg in 1 mL, M = 186.04 g mol⁻¹), degas, and repeat the scan rates. A ferrocene wave of similar height to the analyte wave also indicates a one-electron process.
4. **SWV and LSV.** Record SWV at 20, 40, 60, 80 and 100 Hz, and LSV at 1 and 2 mV s⁻¹, to help count electrons.

## Practical points for the analysis

* **Scan range.** Scan about 500 mV beyond the wave of interest, so that the baseline and the rising part are both well defined for the tangent (onset) construction.
* **Reference choice.** Ferrocene is preferred as the internal reference. Decamethylferrocene can be used when ferrocene overlaps the analyte, but its offset from Fc⁺/Fc must then be known in the same solvent.
* **Frontier orbitals.** Energies come from E_LUMO = −(E_red,onset vs Fc⁺/Fc + 4.8) eV and E_HOMO = −(E_ox,onset vs Fc⁺/Fc + 4.8) eV. Several values for the Fc⁺/Fc level below vacuum are in use (4.8 eV is the most common, 5.1 eV is also used), so state which one you used; `cv-process.sh -E` changes it.
* **Optical gap.** Eg(opt) = 1239.84/λ_onset (nm) from `uv-onset.sh` gives the second frontier orbital when only one redox wave is observed.
