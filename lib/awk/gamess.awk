# =============================================================================
#  gamess.awk - read GAMESS-US and Firefly (PC GAMESS) output files
#
#  usage: awk -v mode=MODE [-v name=NAME] [-v header=1] -f gamess.awk job.log
#  MODE      prints (CSV)
#  summary   name,program,status,opt,charge,mult,natoms,E_final,ZPE,H,G,
#            HOMO_eV,LUMO_eV,gap_eV,dipole_D,n_imag,T_K
#            H and G = E_final + the TOTAL thermal corrections (incl. ZPE)
#  thermo    name,T_K,ZPE_Eh,H_corr_kcal,G_corr_kcal,Cv_cal,S_cal
#  freq      name,mode,freq_cm1,IR_D2_amu_A2   (translations/rotations skipped)
#  charges   name,atom,element,mulliken,lowdin
#  neda      name,component,kcal_mol   (NBO natural energy decomposition:
#            CT, ES, POL, EX or XC, DEF(SE) per fragment, EL, CORE, E)
#  geometry  final Cartesian coordinates (.xyz)
# =============================================================================
BEGIN { H2EV = 27.211386245988; KCAL = 627.509474; if (name == "") name = "mol" }
{ sub(/\r$/, "") }                       # Firefly outputs written on Windows have CRLF endings
/SCF IS UNCONVERGED|SCF DID NOT CONVERGE/ { scfnc = 1 }
/EXECUTION OF GAMESS TERMINATED NORMALLY/  { normal = 1; prog = "GAMESS-US" }
/EXECUTION OF FIREFLY TERMINATED NORMALLY/ { normal = 1; prog = "Firefly" }
/Firefly version|PC GAMESS version/        { prog = "Firefly" }
/GAMESS VERSION =/                         { if (prog == "") prog = "GAMESS-US" }
/EXECUTION OF GAMESS TERMINATED -ABNORMALLY-|ERROR, |\*\*\* ERROR/ { err = 1 }
/EQUILIBRIUM GEOMETRY LOCATED/             { opt = "converged" }
/FAILURE TO LOCATE STATIONARY POINT/       { opt = "not-converged" }
/ MULT  = +[0-9]+ +ICHARG= +-?[0-9]+/ && !cmset { mult = $3; charge = $5; cmset = 1 }
/TOTAL NUMBER OF ATOMS +=/                 { nat = $NF }
/NUMBER OF OCCUPIED ORBITALS \(ALPHA\)/   { nocc = $NF }
/FINAL .*ENERGY IS/                        { for (i = 1; i <= NF; i++) if ($i == "IS") { efinal = $(i + 1); break } }
/^ +TOTAL ENERGY = +-?[0-9]/               { etot = $4 }

# dipole: the line after "DX DY DZ /D/ (DEBYE)"
/DX +DY +DZ +\/D\/ +\(DEBYE\)/ { getline; dip = $4 }

# orbital energies: header line of integers, next line the eigenvalues
/EIGENVECTORS|MOLECULAR ORBITALS$|^ +ALPHA ORBITALS/ { inmo = 1; delete eig; next }
/^ +BETA ORBITALS|END OF (RHF|UHF|ROHF|DFT)|PROPERTY EVALUATION/ { inmo = 0 }
inmo && /^ +[0-9]+( +[0-9]+)* *$/ { nidx = NF; for (k = 1; k <= NF; k++) idx[k] = $k; getline
    if (NF == nidx && $1 ~ /^-?[0-9]+\.[0-9]+$/) for (k = 1; k <= NF; k++) eig[idx[k]] = $k }

# thermochemistry
/THERMOCHEMISTRY AT T= / { T = $4 }
/THE HARMONIC ZERO POINT ENERGY IS/ { getline; zpe = $1 }
/KCAL\/MOL +KCAL\/MOL +KCAL\/MOL CAL\/MOL-K/ { inkc = 1; next }
inkc && /^ +TOTAL / { th = $3; tg = $4; tcv = $5; ts = $7; inkc = 0; nT++; tT[nT] = T; tZ[nT] = zpe; tH[nT] = th; tG[nT] = tg; tC[nT] = tcv; tS[nT] = ts }

# frequencies (blocks of five modes)
/^ +FREQUENCY: / { nb = 0; for (k = 2; k <= NF; k++) { if ($k ~ /^I$/) { fr[nf] = -fr[nf]; continue } nf++; fr[nf] = $k; nb++ } fstart = nf - nb }
/^ +IR INTENSITY: / { for (k = 3; k <= NF; k++) ir[fstart + k - 2] = $k }

# charges (last Mulliken/Lowdin table)
/TOTAL MULLIKEN AND LOWDIN ATOMIC POPULATIONS/ { inq = 1; nq = 0; getline; next }
inq && /^ +[0-9]+ +[A-Z][A-Za-z]* +[0-9.]+ +-?[0-9.]+ +[0-9.]+ +-?[0-9.]+/ { nq++; qe[nq] = $2; qm[nq] = $4; ql[nq] = $6; next }
inq && nq && !/^ +[0-9]+ / { inq = 0 }

# geometry: input orientation in bohr (single points) and every optimisation step in angstrom
function elname(e) { return toupper(substr(e, 1, 1)) tolower(substr(e, 2)) }
/ATOM +ATOMIC +COORDINATES \(BOHR\)/ && !angs { inb = 1; nb_ = 0; getline; next }
inb { if (NF != 5) { inb = 0; ng = nb_; next } nb_++; ge[nb_] = elname($1); gx[nb_] = $3 * 0.529177210903; gy[nb_] = $4 * 0.529177210903; gz[nb_] = $5 * 0.529177210903; next }
/COORDINATES OF ALL ATOMS ARE \(ANGS\)/ { angs = 1 }
/COORDINATES OF ALL ATOMS ARE \(ANGS\)/ { ing = 1; ng = 0; getline; getline; next }
ing { if (NF != 5) { ing = 0; next } ng++; ge[ng] = elname($1); gx[ng] = $3; gy[ng] = $4; gz[ng] = $5; next }

# NBO natural energy decomposition analysis (summary block)
/Natural Energy Decomposition Analysis \(Summary\)/ { inneda = 1; nn = 0; ndef = 0; next }
inneda && / CT = / { nn++; nk[nn] = "CT"; nv[nn] = $NF }
inneda && / ES = / { nn++; nk[nn] = "ES"; nv[nn] = $NF }
inneda && / POL = / { nn++; nk[nn] = "POL"; nv[nn] = $NF }
inneda && / (EX|XC) = / { for (k = 1; k <= NF; k++) if ($k == "EX" || $k == "XC") { nn++; nk[nn] = $k; nv[nn] = $NF } }
inneda && /DEF\(SE\) = / { line = $0; sub(/.*DEF\(SE\) = */, "", line); gsub(/[()]/, " ", line); split(line, w, " ")
    ndef++; nn++; nk[nn] = "DEF_frag" ndef; nv[nn] = w[1]; nn++; nk[nn] = "SE_frag" ndef; nv[nn] = w[2] }
inneda && /^ +E = / { nn++; nk[nn] = "E_sum"; nv[nn] = $NF }
inneda && /Electrical .*:/ { nn++; nk[nn] = "EL(ES+POL+SE)"; nv[nn] = $NF }
inneda && /Charge Transfer .*:/ { nn++; nk[nn] = "CT_total"; nv[nn] = $NF }
inneda && /Core .*:/ { c = $0; sub(/ *:.*/, "", c); gsub(/ /, "", c); nn++; nk[nn] = c; nv[nn] = $NF }
inneda && /Total Interaction .*:/ { nn++; nk[nn] = "Total_interaction"; nv[nn] = $NF; inneda = 0 }

END {
    status = (err ? "error" : (normal ? "normal" : "incomplete")); if (scfnc) status = status "/SCF-not-converged"
    if (mode == "summary") {
        e = (efinal != "" ? efinal : etot); nim = 0
        for (i = 1; i <= nf; i++) if (fr[i] < 0 && fr[i] < -20) nim++
        homo = (nocc in eig ? eig[nocc] * H2EV : ""); lumo = ((nocc + 1) in eig ? eig[nocc + 1] * H2EV : "")
        H = (th != "" ? sprintf("%.6f", e + th / KCAL) : ""); G = (tg != "" ? sprintf("%.6f", e + tg / KCAL) : "")
        if (header) print "name,program,status,opt,charge,mult,natoms,E_final,ZPE,H,G,HOMO_eV,LUMO_eV,gap_eV,dipole_D,n_imag,T_K"
        printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n", name, prog, status, opt, charge, mult, (nat ? nat : ng), e, zpe, H, G,
            (homo != "" ? sprintf("%.4f", homo) : ""), (lumo != "" ? sprintf("%.4f", lumo) : ""),
            (homo != "" && lumo != "" ? sprintf("%.4f", lumo - homo) : ""), dip, (nf ? nim : ""), T
    }
    else if (mode == "thermo") {
        if (header) print "name,T_K,ZPE_Eh,H_corr_kcal,G_corr_kcal,Cv_cal_molK,S_cal_molK"
        for (i = 1; i <= nT; i++) print name "," tT[i] "," tZ[i] "," tH[i] "," tG[i] "," tC[i] "," tS[i]
    }
    else if (mode == "freq") {
        if (header) print "name,mode,freq_cm1,IR_D2_amu_A2"
        k = 0; for (i = 1; i <= nf; i++) if (fr[i] > 20 || fr[i] < -20) print name "," ++k "," fr[i] "," ir[i]
    }
    else if (mode == "charges") {
        if (header) print "name,atom,element,mulliken,lowdin"
        for (i = 1; i <= nq; i++) print name "," i "," qe[i] "," qm[i] "," ql[i]
    }
    else if (mode == "neda") {
        if (header) print "name,component,kcal_mol"
        for (i = 1; i <= nn; i++) print name "," nk[i] "," nv[i]
    }
    else if (mode == "geometry") {
        print ng; print name " final geometry (" status ")"
        for (i = 1; i <= ng; i++) printf "%-3s %14.8f %14.8f %14.8f\n", ge[i], gx[i], gy[i], gz[i]
    }
    else { print "unknown mode " mode > "/dev/stderr"; exit 1 }
}
