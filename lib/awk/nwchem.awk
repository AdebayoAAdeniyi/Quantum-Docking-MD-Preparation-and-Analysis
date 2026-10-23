# =============================================================================
#  nwchem.awk - read NWChem 6.x/7.x output
#  usage: awk -v mode=MODE [-v name=NAME] [-v header=1] -f nwchem.awk job.log
#  summary   name,status,opt_converged,charge,mult,E_final,ZPE,H,G,HOMO_eV,LUMO_eV,gap_eV,n_imag,lowest_freq,T_K
#            (G = E + H_corr - T*S, since NWChem prints H and S but not G; 298.15 K if computed)
#  thermo    name,T_K,ZPE,E_corr,H_corr,S_cal_molK,Cv_cal_molK,G_corr   (au unless stated)
#  freq      name,mode,freq_cm1,IR_km_mol     (|freq| < 20 cm-1 skipped)
#  neb       name,bead,E_Eh,dE_kcal_mol        (last "Path Energy" table)
#  bsse      name,E_super,BSSE_error,E_corrected
#  geometry  last "Output coordinates in angstroms" block as .xyz
# =============================================================================
BEGIN { H2EV = 27.211386245988; KCAL = 627.509474; if (name == "") name = "mol" }
/Total times +cpu:/ { normal = 1 }
/There is an error|MPI_ABORT|ARMCI Error/ { err = 1 }
/Optimization converged/ { opt = "yes" }
/Optimization failed|maximum number of steps/ { if (opt == "") opt = "no" }
/^ +task +dft +(optimize|saddle)/ { if (opt == "") opt = "no" }
/^ +charge +=/ && chg == "" { chg = $3 + 0 }
/^ +mult +[0-9]+/ { mult = $2 }
/Spin multiplicity:/ { mult = $NF }
/Total DFT energy =|Total SCF energy =/ { e = $NF }
/^ +HOMO +=/ { homo = $3 } /^ +LUMO +=/ { lumo = $3 }

# thermochemistry, one block per temperature
/^ Temperature +=/ { T = $3; sub(/K/, "", T) }
/Zero-Point correction to Energy/ { z = $(NF - 1) }
/Thermal correction to Energy/ { ec = $(NF - 1) }
/Thermal correction to Enthalpy/ { hc = $(NF - 1) }
/^ Total Entropy +=/ { s = $4 }
/Cv \(constant volume heat capacity\) =/ { cv = $7; nT++; tT[nT] = T; tZ[nT] = z; tE[nT] = ec; tH[nT] = hc; tS[nT] = s; tC[nT] = cv }

# projected frequencies and IR intensities (last table)
/Projected Infra Red Intensities/ { inir = 1; nf = 0; getline; getline; next }
inir && /^ +[0-9]+ +-?[0-9.]+ +\|\|/ { nf++; fr[nf] = $2; ir[nf] = $6; next }
inir && /^ +-{10,}/ { inir = 0 }

# NEB path energies (last table)
/neb: Path Energy #/ { innb = 1; nb = 0; next }
innb && /^ neb: +[0-9]+ +-[0-9.]+ *$/ { nb++; bi[nb] = $2; be[nb] = $3; next }
innb && nb && !/^ neb: +[0-9]+ / { innb = 0 }

# counterpoise
/BSSE error =/ { bsse = $NF } /Supermolecular energy =/ { esup = $NF } /Corrected energy =/ { ecor = $NF }

# geometry
/Output coordinates in angstroms/ { ing = 1; ng = 0; skip = 3; next }
ing { if (skip-- > 0) next; if (NF < 6) { ing = 0; next }
      ng++; el = $2; sub(/^[Bb][Qq]/, "", el); ge[ng] = toupper(substr(el,1,1)) tolower(substr(el,2)); gx[ng] = $4; gy[ng] = $5; gz[ng] = $6; next }

END {
    st = (err ? "error" : (normal ? "normal" : "incomplete"))
    if (mode == "summary") {
        nim = 0; low = ""
        for (i = 1; i <= nf; i++) if (fr[i] > 20 || fr[i] < -20) { if (fr[i] < 0) nim++; if (low == "" || fr[i] < low) low = fr[i] }
        # report 298.15 K when it was computed, otherwise the last temperature
        for (i = 1; i <= nT; i++) if (tT[i] + 0 == 298.15) { T = tT[i]; z = tZ[i]; hc = tH[i]; s = tS[i] }
        H = (hc != "" ? sprintf("%.6f", e + hc) : ""); G = (hc != "" ? sprintf("%.6f", e + hc - T * s / 1000 / KCAL) : "")
        if (header) print "name,status,opt_converged,charge,mult,E_final,ZPE,H,G,HOMO_eV,LUMO_eV,gap_eV,n_imag,lowest_freq,T_K"
        printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n", name, st, opt, chg, mult, e, z, H, G,
            (homo != "" ? sprintf("%.4f", homo * H2EV) : ""), (lumo != "" ? sprintf("%.4f", lumo * H2EV) : ""),
            (homo != "" && lumo != "" ? sprintf("%.4f", (lumo - homo) * H2EV) : ""), (nf ? nim : ""), low, T
    } else if (mode == "thermo") {
        if (header) print "name,T_K,ZPE,E_corr,H_corr,S_cal_molK,Cv_cal_molK,G_corr"
        for (i = 1; i <= nT; i++) printf "%s,%s,%s,%s,%s,%s,%s,%.6f\n", name, tT[i], tZ[i], tE[i], tH[i], tS[i], tC[i], tH[i] - tT[i] * tS[i] / 1000 / KCAL
    } else if (mode == "freq") {
        if (header) print "name,mode,freq_cm1,IR_km_mol"
        k = 0; for (i = 1; i <= nf; i++) if (fr[i] > 20 || fr[i] < -20) print name "," ++k "," fr[i] "," ir[i]
    } else if (mode == "neb") {
        if (header) print "name,bead,E_Eh,dE_kcal_mol"
        for (i = 1; i <= nb; i++) printf "%s,%s,%s,%.2f\n", name, bi[i], be[i], (be[i] - be[1]) * KCAL
    } else if (mode == "bsse") {
        if (header) print "name,E_super,BSSE_error,E_corrected,BSSE_kcal_mol"
        printf "%s,%s,%s,%s,%.3f\n", name, esup, bsse, ecor, bsse * KCAL
    } else if (mode == "geometry") {
        print ng; print name " final geometry (" st ")"
        for (i = 1; i <= ng; i++) printf "%-3s %14.8f %14.8f %14.8f\n", ge[i], gx[i], gy[i], gz[i]
    } else { print "unknown mode " mode > "/dev/stderr"; exit 1 }
}
