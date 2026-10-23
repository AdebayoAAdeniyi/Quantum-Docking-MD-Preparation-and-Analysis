# =============================================================================
#  orca.awk - read ORCA 4.2 / 5 / 6 output files
#
#  usage: awk -v mode=MODE [-v name=NAME] [-v header=1] -f orca.awk job.out
#
#  MODE       prints (CSV)
#  summary    name,status,opt_converged,charge,mult,natoms,E_final,ZPE,H,G,
#             HOMO_eV,LUMO_eV,gap_eV,dipole_D,n_imag,lowest_freq,T_K
#  thermo     name,T_K,E_el,ZPE,thermal_corr,H,TS,G  (one row per temperature)
#  charges    name,atom,element,mulliken,loewdin,hirshfeld,chelpg,spin_mulliken
#  freq       name,mode,freq_cm1,IR_km_mol,Raman_activity
#  uv         name,state,multiplicity,E_eV,lambda_nm,f,main_transitions
#             (TD-DFT/CIS; ORCA 5 and ORCA 6 table layouts)
#  socuv      name,state,E_eV,lambda_nm,f   (spin-orbit corrected spectrum)
#  socme      name,T,S,SOCME_cm1  |<T|HSO|S>| = sqrt(sum of Re^2 + Im^2 over x,y,z)
#  casscf     name,state,root,mult,dE_CASSCF_eV,dE_NEVPT2_eV  + active occupations
#  neb        name,image,E_Eh,dE_kcal_mol,max_Fp,marker  (final path summary)
#  esd        name,kind,rate_s-1,FC_percent,HT_percent
#  nmr        name,atom,element,shielding_iso_ppm,anisotropy_ppm
#  geometry   final Cartesian coordinates as an .xyz file
#  Atoms are numbered from 1 here (ORCA counts from 0).
# =============================================================================
function trim(s) { gsub(/^[ \t]+|[ \t]+$/, "", s); return s }
BEGIN { H2EV = 27.211386245988; if (name == "") name = "mol" }

/ORCA TERMINATED NORMALLY/                 { normal = 1 }
/error termination|[Aa]borting the run|ORCA finished by error/ { err = 1 }
/THE OPTIMIZATION HAS CONVERGED/            { optconv = "yes" }
/did not converge but reached the maximum number/ { optconv = "no" }
/Geometry Optimization Run/                 { if (optconv == "") optconv = "no" }
/Total Charge +Charge +\.\.\.\./            { charge = $NF }
/Multiplicity +Mult +\.\.\.\./              { mult = $NF }
/Number of atoms +\.\.\./                   { nat = $NF }
/^FINAL SINGLE POINT ENERGY/                { efinal = $5 }
/SCF NOT CONVERGED/                         { scfnc = 1 }
/Magnitude \(Debye\) +:/                    { dip = $NF }

# ---------- geometry ----------
/^CARTESIAN COORDINATES \(ANGSTROEM\)/ { ingeo = 1; ng = 0; getline; next }
ingeo { if (NF < 4) { ingeo = 0; next } ng++; ge[ng] = $1; gx[ng] = $2; gy[ng] = $3; gz[ng] = $4; next }

# ---------- orbitals (alpha / closed shell) ----------
/^ORBITAL ENERGIES/ { inorb = 1; homo = ""; lumo = ""; next }
inorb && /SPIN DOWN ORBITALS/ { inorb = 0 }
inorb && /^ +[0-9]+ +[0-9.]+ +-?[0-9.]+ +-?[0-9.]+/ { if ($2 > 0) homo = $4; else if (lumo == "") lumo = $4 }
inorb && /^-{5,}|^\*/ && homo != "" { inorb = 0 }

# ---------- thermochemistry (one block per temperature) ----------
/^Temperature +\.\.\./              { T = $3 }
/^Electronic energy +\.\.\./        { eel = $4 }
/^Zero point energy +\.\.\./        { zpe = $5 }
/^Total thermal correction/         { tcor = $4 }
/^Total Enthalpy +\.\.\./           { hh = $4 }
/^Final entropy term +\.\.\./       { ts = $5 }
/^Final Gibbs free energy +\.\.\./  { gg = $6; nT++; tT[nT] = T; tE[nT] = eel; tZ[nT] = zpe; tC[nT] = tcor; tH[nT] = hh; tS[nT] = ts; tG[nT] = gg }

# ---------- vibrations ----------
/^VIBRATIONAL FREQUENCIES/ { invib = 1; nf = 0; next }
invib && /^ +[0-9]+: +-?[0-9.]+ cm\*\*-1/ { m = $1 + 0; fr[m] = $2; if (m > nfm) nfm = m; next }
invib && /^NORMAL MODES/ { invib = 0 }
/^IR SPECTRUM/    { inir = 1; nir = 0; next }
inir && /^ +[0-9]+: +[0-9.-]+ +[0-9.]+ +[0-9.]+/ { ir[$1 + 0] = $4; nir++; next }
inir && nir && /^ *$/ { inir = 0 }
/^RAMAN SPECTRUM/ { inra = 1; next }
inra && /^ +[0-9]+: +[0-9.-]+ +[0-9.]+/ { ra[$1 + 0] = $3; next }

# ---------- excited states ----------
/^STATE +[0-9]+: +E= / { cur = sprintf("%.0f", $8); ntr[cur] = 0; instate = 1; next }
instate && /^ +[0-9]+[ab] +-> +[0-9]+[ab] +: +[0-9.]+/ {
    if ($5 >= 0.10) tr[cur] = tr[cur] (tr[cur] == "" ? "" : "; ") $1 "->" $3 " (" sprintf("%.0f", 100 * $5) "%)"
    next }
instate && /^$/ { instate = 0 }
/^ +ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS/ && !soc && !absdone { inabs = 1; next }
/SPIN ORBIT CORRECTED ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC/ { insoc = 1; next }
inabs && /^ +[0-9]+-[0-9]*[A-Za-z]* +-> +[0-9]+-[0-9]+/ {          # ORCA 6
    split($3, a, "-"); m = a[2]; gsub(/[^0-9]/, "", m)
    ns++; st[ns] = a[1]; sm[ns] = (m == 3 ? "Triplet" : m == 1 ? "Singlet" : m); se[ns] = $4; sl[ns] = $6; sf[ns] = $7
    sk[ns] = sprintf("%.0f", $5); next }
inabs && /^ +[0-9]+ +[0-9.]+ +[0-9.]+ +[0-9.]+ +[0-9.]+/ {          # ORCA 4/5
    ns++; st[ns] = $1; sm[ns] = "Singlet"; se[ns] = sprintf("%.4f", $2 / 8065.544); sl[ns] = $3; sf[ns] = $4
    sk[ns] = sprintf("%.0f", $2); next }
inabs && /spin forbidden \(mult=3\)/ {
    ns++; st[ns] = $1; sm[ns] = "Triplet"; se[ns] = sprintf("%.4f", $2 / 8065.544); sl[ns] = $3; sf[ns] = 0; sk[ns] = sprintf("%.0f", $2); next }
inabs && ns && /^$/ { inabs = 0; absdone = 1 }
insoc && /^ +[0-9]+ +[0-9]+ +[0-9.]+ +[0-9.]+ +[0-9.]+/ { nsoc++; so_s[nsoc] = $2; so_e[nsoc] = sprintf("%.4f", $3 / 8065.544); so_l[nsoc] = $4; so_f[nsoc] = $5; next }
insoc && /^ +[0-9]+-[0-9.A-Za-z]* +-> +[0-9]+/ { nsoc++; so_s[nsoc] = $3; so_e[nsoc] = $4; so_l[nsoc] = $6; so_f[nsoc] = $7; next }
insoc && nsoc && /^$/ { insoc = 0 }

# ---------- spin-orbit coupling matrix elements ----------
/CALCULATED SOCME BETWEEN TRIPLETS AND SINGLETS/ { insocme = 1; nsm = 0; next }
insocme && /^ +[0-9]+ +[0-9]+ +\(/ {
    line = $0; gsub(/[(),]/, " ", line); n = split(line, w, " ")
    nsm++; smT[nsm] = w[1]; smS[nsm] = w[2]; s2 = 0; for (k = 3; k <= 8; k++) s2 += w[k] * w[k]; smV[nsm] = sqrt(s2); next }
insocme && nsm && /^ *-*$/ { insocme = 0 }

# ---------- CASSCF / NEVPT2 ----------
/^ +N\(occ\)= / { occ = $0; sub(/^ +N\(occ\)= +/, "", occ) }
/^Final CASSCF energy/ { ecas = $5 }
/^SA-CASSCF TRANSITION ENERGIES/ { incas = 1; crow = 0; next }
/NEVPT2 TRANSITION ENERGIES/ && !/CORRECTION/ { innev = 1; crow = 0; next }
(incas || innev) && /^ +[0-9]+: +[0-9]+ +[0-9]+ +-?[0-9.]+ +-?[0-9.]+/ {
    k = $1 + 0; crow++; if (incas) { cR[k] = $2; cM[k] = $3; cE[k] = $5; if (k > ncas) ncas = k } else nE[k] = $5; next }
(incas || innev) && crow && /^ *$/ { incas = 0; innev = 0 }

# ---------- NEB path summary (the last table printed) ----------
/^ *Image +E\(Eh\) +dE\(kcal\/mol\)/ { innb = 1; nnb = 0; next }
innb && /^ +([0-9]+|TS) +-?[0-9]+\.[0-9]+ +-?[0-9.]+/ {
    nnb++; nbI[nnb] = $1; nbE[nnb] = $2; nbD[nnb] = $3; nbF[nnb] = $4; nbM[nnb] = ($0 ~ /<= CI/ ? "CI" : $0 ~ /<= TS/ ? "TS" : ""); next }
innb && nnb && !/^ +([0-9]+|TS) / { innb = 0 }

# ---------- excited state dynamics (ESD) rates ----------
/The calculated (fluorescence|phosphorescence|internal conversion|intersystem crossing) rate constant is/ {
    ne++; kind = $0; sub(/.*The calculated /, "", kind); sub(/ rate constant.*/, "", kind); esK[ne] = kind
    for (k = 1; k <= NF; k++) if ($k ~ /^[0-9.]+[eE][-+]?[0-9]+$/) esR[ne] = $k
    getline; if ($0 ~ /from FC/) { split($0, q, "%"); n1 = split(q[1], u, " "); esF[ne] = u[n1]; n2 = split(q[2], v, " "); esH[ne] = v[n2] } }

# ---------- charges ----------
/^MULLIKEN ATOMIC CHARGES/ { inmu = 1; nmu = 0; getline; next }
inmu { if ($0 !~ /^ +[0-9]+ +[A-Za-z]+ *:/) { inmu = 0; next } nmu++; el[nmu] = $2; mu[nmu] = $4; if (NF >= 5) msp[nmu] = $5; next }
/^LOEWDIN ATOMIC CHARGES/ { inlo = 1; nlo = 0; getline; next }
inlo { if ($0 !~ /^ +[0-9]+ +[A-Za-z]+ *:/) { inlo = 0; next } lo[++nlo] = $4; next }
/^HIRSHFELD ANALYSIS/ { inhi = 1; nhi = 0; next }
inhi && /^ +[0-9]+ +[A-Za-z]+ +-?[0-9.]+/ { hi[++nhi] = $3; next }
inhi && /^ +TOTAL/ { inhi = 0 }
/CHELPG Charges/ { inch = 1; nch = 0; getline; next }
inch { if ($0 !~ /^ +[0-9]+ +[A-Za-z]+ *:/) { inch = 0; next } ch[++nch] = $NF; next }

# ---------- NMR ----------
/^CHEMICAL SHIELDING SUMMARY/ { innmr = 1; nn = 0; next }
innmr && /^ +[0-9]+ +[A-Z][a-z]? +-?[0-9.]+ +-?[0-9.]+/ { nn++; na[nn] = $1 + 1; ne_[nn] = $2; ni[nn] = $3; nan[nn] = $4; next }
innmr && nn && /^$/ { innmr = 0 }

END {
    status = (err ? "error" : (normal ? "normal" : "incomplete")); if (scfnc) status = status "/SCF-not-converged"
    if (mode == "summary") {
        nim = 0; low = ""
        for (i = 0; i <= nfm; i++) if (i in fr && fr[i] != 0) { if (fr[i] < 0) nim++; if (low == "" || fr[i] < low) low = fr[i] }
        gap = (homo != "" && lumo != "" ? sprintf("%.4f", lumo - homo) : "")
        if (header) print "name,status,opt_converged,charge,mult,natoms,E_final,ZPE,H,G,HOMO_eV,LUMO_eV,gap_eV,dipole_D,n_imag,lowest_freq,T_K"
        print name "," status "," optconv "," charge "," mult "," (nat ? nat : ng) "," efinal "," zpe "," hh "," gg "," homo "," lumo "," gap "," dip "," (nfm ? nim : "") "," low "," T
    }
    else if (mode == "thermo") {
        if (header) print "name,T_K,E_el,ZPE,thermal_corr,H,TS,G"
        for (i = 1; i <= nT; i++) print name "," tT[i] "," tE[i] "," tZ[i] "," tC[i] "," tH[i] "," tS[i] "," tG[i]
    }
    else if (mode == "charges") {
        if (header) print "name,atom,element,mulliken,loewdin,hirshfeld,chelpg,spin_mulliken"
        for (i = 1; i <= nmu; i++) print name "," i "," el[i] "," mu[i] "," lo[i] "," hi[i] "," ch[i] "," msp[i]
    }
    else if (mode == "freq") {
        if (header) print "name,mode,freq_cm1,IR_km_mol,Raman_activity"
        for (i = 0; i <= nfm; i++) if (i in fr && fr[i] != 0) print name "," i "," fr[i] "," ir[i] "," ra[i]
    }
    else if (mode == "uv") {
        if (header) print "name,state,multiplicity,E_eV,lambda_nm,f,main_transitions(>=10%)"
        for (i = 1; i <= ns; i++) print name "," st[i] "," sm[i] "," se[i] "," sl[i] "," sf[i] ",\"" tr[sk[i]] "\""
    }
    else if (mode == "socuv") {
        if (header) print "name,state,E_eV,lambda_nm,f"
        for (i = 1; i <= nsoc; i++) print name "," so_s[i] "," so_e[i] "," so_l[i] "," so_f[i]
    }
    else if (mode == "socme") {
        if (header) print "name,T_root,S_root,SOCME_cm1"
        for (i = 1; i <= nsm; i++) printf "%s,%s,%s,%.2f\n", name, smT[i], smS[i], smV[i]
    }
    else if (mode == "casscf") {
        if (header) print "name,state,root,mult,dE_CASSCF_eV,dE_NEVPT2_eV,E_CASSCF_Eh,active_occupations"
        for (i = 1; i <= ncas; i++) print name "," i "," cR[i] "," cM[i] "," cE[i] "," nE[i] "," ecas ",\"" occ "\""
        n = split(occ, o, " "); for (i = 1; i <= n; i++) if (o[i] < 0.02 || o[i] > 1.98)
            printf "warning: %s active orbital %d has occupation %s; it is barely correlated, consider removing it from the active space\n", name, i, o[i] > "/dev/stderr"
    }
    else if (mode == "neb") {
        if (header) print "name,image,E_Eh,dE_kcal_mol,max_Fp,marker"
        for (i = 1; i <= nnb; i++) print name "," nbI[i] "," nbE[i] "," nbD[i] "," nbF[i] "," nbM[i]
    }
    else if (mode == "esd") {
        if (header) print "name,kind,rate_s-1,FC_percent,HT_percent"
        for (i = 1; i <= ne; i++) print name "," esK[i] "," esR[i] "," esF[i] "," esH[i]
    }
    else if (mode == "nmr") {
        if (header) print "name,atom,element,shielding_iso_ppm,anisotropy_ppm"
        for (i = 1; i <= nn; i++) print name "," na[i] "," ne_[i] "," ni[i] "," nan[i]
    }
    else if (mode == "geometry") {
        print ng; print name " final geometry (" status ")"
        for (i = 1; i <= ng; i++) printf "%-3s %14.8f %14.8f %14.8f\n", ge[i], gx[i], gy[i], gz[i]
    }
    else { print "unknown mode: " mode > "/dev/stderr"; exit 1 }
}
