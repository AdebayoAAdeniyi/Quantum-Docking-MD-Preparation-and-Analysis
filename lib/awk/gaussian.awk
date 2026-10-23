# =============================================================================
#  gaussian.awk - read Gaussian 09/16 output (.log/.out) files
#
#  usage: awk -v mode=MODE [-v name=NAME] -f gaussian.awk file.log
#
#  MODE        prints (comma separated, one header line when -v header=1)
#  summary     name,status,charge,mult,natoms,E_scf,E_post,post_method,
#              ZPE,E_thermal,H,G,Cv,S,
#              HOMO_eV,LUMO_eV,gap_eV,dipole_D,n_imag,lowest_freq,T_K
#              (HOMO/LUMO are alpha orbitals; for open shells read the log)
#  thermo      name,T_K,P_atm,E_scf,ZPE_corr,H_corr,G_corr,E0,H,G,Cv,S
#              (one row per temperature, for freq(...) Temperature=... jobs
#               and for logs that contain several frequency steps)
#  charges     name,atom,element,mulliken,apt,npa,hirshfeld,cm5,spin
#  uv          name,state,multiplicity,E_eV,lambda_nm,f,main_transitions
#  freq        name,mode,freq_cm1,IR_km_mol,Raman_A4_amu
#  nmr         name,atom,element,shielding_iso,anisotropy
#  polar       name,alpha_iso_au,alpha_aniso_au,beta_x,beta_y,beta_z,
#              beta_tot_au,beta_tot_1e-30esu
#  hyperfine   name,atom,isotope,fermi_MHz,fermi_Gauss
#  jcoupling   name,atom_i,el_i,atom_j,el_j,FC,SD,PSO,DSO,J_total (Hz)
#              [-v pairs="1-2,3-7"] [-v jmin=0.5]
#  geometry    final Cartesian coordinates as an .xyz file
#
#  "Last value wins": for multi-step (--Link1--) logs split the file first
#  with split-link1-log.sh if you need every step separately.
#  Energies are in Hartree unless the column name says otherwise.
# =============================================================================
function trim(s) { gsub(/^[ \t]+|[ \t]+$/, "", s); return s }
function dnum(s) { gsub(/D/, "E", s); return sprintf("%.12g", s + 0) }   # Fortran D exponents, full precision
BEGIN {
    H2EV = 27.211386245988
    BAU2ESU = 8.6393e-3       # 1 a.u. of beta = 8.6393e-33 esu = 8.6393e-3 x 1e-30 esu
    if (name == "") name = "mol"
    nT = 0
}

# ---------------- general ----------------
/Normal termination of Gaussian/ { normal++ }
/Error termination/             { error++ }
/ Charge = .* Multiplicity = / && !cm_done { charge = $3; mult = $6; cm_done = 1 }
/^ NAtoms=/ { nat = $2 }
/SCF Done:/ { escf = $5 }
# correlated energies: the highest level printed is reported as E_post
/EUMP2 =/                  { split($0, a, "EUMP2 =");  epost = dnum(trim(a[2])); lpost = "MP2" }
/UMP4\(SDTQ\)=/            { split($0, a, "UMP4\\(SDTQ\\)="); epost = dnum(trim(a[2])); lpost = "MP4" }
/^ CCSD\(T\)= /             { epost = dnum($2); lpost = "CCSD(T)" }

# ---------------- geometry (last orientation block) ----------------
/(Standard|Input) orientation:/ { ingeo = 1; ngeo = 0; skip = 4; next }
ingeo {
    if (skip > 0) { skip--; next }
    if ($0 ~ /-----/) { ingeo = 0; nat = ngeo; next }
    ngeo++; gz[ngeo] = $2; gx[ngeo] = $4; gy[ngeo] = $5; gzz[ngeo] = $6; next
}

# ---------------- orbital energies ----------------
/Alpha  occ\. eigenvalues/ {
    if (lastblock != "occ") { homo = ""; lumo = "" }      # a new population block
    homo = $NF; lastblock = "occ"
}
/Alpha virt\. eigenvalues/ {
    if (lastblock == "occ") lumo = $5
    lastblock = "virt"
}

# ---------------- dipole ----------------
/ Dipole moment \(field-independent basis, Debye\):/ { getline; dip = $8 }

# ---------------- thermochemistry ----------------
/^ Temperature +[0-9.]+ Kelvin\.  Pressure/ { T = $2; P = $5 }
/Zero-point correction=/                       { zpe = $3 }
/Thermal correction to Enthalpy=/              { hcorr = $5 }
/Thermal correction to Gibbs Free Energy=/     { gcorr = $7 }
/Sum of electronic and zero-point Energies=/   { e0 = $7 }
/Sum of electronic and thermal Energies=/      { eth = $7 }
/Sum of electronic and thermal Enthalpies=/    { hh = $7 }
/Sum of electronic and thermal Free Energies=/ { gg = $8 }
/E \(Thermal\) +CV +S/ { getline; getline; cv = $3; ss = $4
    nT++; tT[nT] = T; tP[nT] = P; tZ[nT] = zpe; tH[nT] = hcorr; tG[nT] = gcorr
    tE0[nT] = e0; tHH[nT] = hh; tGG[nT] = gg; tCv[nT] = cv; tS[nT] = ss; tE[nT] = escf }

# ---------------- frequencies ----------------
/^ +Frequencies -- / { if (!freqstarted) { nf = 0; freqstarted = 1 }
                       for (i = 3; i <= NF; i++) fr[++nf] = $i; nfi = nf - (NF - 2); nfr = nfi }
/^ +IR Inten +-- /   { for (i = 4; i <= NF; i++) ir[++nfi] = $i }
/^ +Raman Activ +-- / { for (i = 4; i <= NF; i++) ra[++nfr] = $i }
/Harmonic frequencies \(cm\*\*-1\)/ { freqstarted = 0 }     # a new frequency job restarts the list

# ---------------- charges ----------------
/^ (Mulliken charges|Mulliken charges and spin densities|Mulliken atomic charges):/ { inmul = 1; nm = 0; getline; next }
inmul { if ($0 ~ /Sum of Mulliken/) { inmul = 0; next }
        nm++; el[nm] = $2; mul[nm] = $3; if (NF >= 4) spin[nm] = $4; next }
/^ APT charges:/ { inapt = 1; na = 0; getline; next }
inapt { if ($0 ~ /Sum of APT/) { inapt = 0; next } apt[++na] = $3; next }
# NBO prints the total-density summary first, then alpha and beta spin
# summaries for open shells: keep only the first summary of each NBO run
/N A T U R A L   A T O M I C   O R B I T A L/ { npa_lock = 0 }
/Summary of Natural Population Analysis:/ && !npa_lock { innpa = 1; npa_lock = 1; nn = 0; skipn = 5; next }
innpa { if (skipn > 0) { skipn--; next }
        if ($0 ~ /=====/) { innpa = 0; next }
        npa[++nn] = $3; if (el[nn] == "") el[nn] = $1; next }
/Hirshfeld charges, spin densities, dipoles, and CM5 charges|Hirshfeld spin densities, charges and dipoles/ {
        inhir = 1; nh = 0; getline; next }
inhir { if ($1 == "Tot" || $0 ~ /^ *$/ || $0 ~ /Hirshfeld charges with hydrogens/) { inhir = 0; next }
        nh++; hir[nh] = $3; hspin[nh] = $4; if (NF >= 8) cm5[nh] = $8; if (el[nh] == "") el[nh] = $2; next }

# ---------------- excited states ----------------
/^ Excited State +[0-9]+:/ {
    ns++; sub(/:/, "", $3); st[ns] = $3
    split($4, m, "-"); smult[ns] = m[1]; sev[ns] = $5; snm[ns] = $7
    f = $9; sub(/f=/, "", f); sf[ns] = f; strans[ns] = ""; nbest[ns] = 0; next
}
ns && /^ +[0-9]+[AB]? *-> *[0-9]+[AB]? +-?[0-9.]+$/ {
    line = $0; gsub(/->/, " -> ", line); n2 = split(line, w, " ")
    c = w[n2]; pct = 200 * c * c                       # 2c^2 for closed shell, as %
    if (pct >= 10) strans[ns] = strans[ns] (strans[ns] == "" ? "" : "; ") w[1] "->" w[3] " (" sprintf("%.0f", pct) "%)"
}

# ---------------- NMR ----------------
/Isotropic = .*Anisotropy =/ { nmr_n++; nmr_at[nmr_n] = $1; nmr_el[nmr_n] = $2; nmr_iso[nmr_n] = $5; nmr_an[nmr_n] = $8 }
/Magnetic shielding tensor/ { nmr_n = 0 }

# ---------------- polarisability ----------------
/^ +iso +-?[0-9.]+D[+-][0-9]+/   && inalpha { a_iso = dnum($2) }
/^ +aniso +-?[0-9.]+D[+-][0-9]+/ && inalpha { a_aniso = dnum($2) }
/Dipole polarizability, Alpha/ { inalpha = 1; inbeta = 0 }
/First dipole hyperpolarizability, Beta/ { inalpha = 0; if (!beta_done) inbeta = 1 }   # first block = input orientation
inbeta && /^ +\|\| +-?[0-9.]+D/ { bpar = dnum($2) }
/Isotropic polarizability for W= +0\.0+ / && a_iso == "" { a_iso_simple = $6 }
inbeta && /^ +(xxx|xxy|yxy|yyy|xxz|yxz|yyz|zxz|zyz|zzz) +-?[0-9.]+D/ { B[$1] = dnum($2); nB++; if ($1 == "zzz") { inbeta = 0; beta_done = 1 } }

# ---------------- spin-spin couplings (NMR=SpinSpin), lower triangle blocks ----------------
/contribution to J \(Hz\):|Total nuclear spin-spin coupling J \(Hz\):/ {
    jt = ($0 ~ /Fermi/ ? "FC" : $0 ~ /Spin-dipolar/ ? "SD" : $0 ~ /Paramagnetic/ ? "PSO" : $0 ~ /Diamagnetic/ ? "DSO" : "J")
    inj = 1; next }
inj {
    if ($0 ~ /^ +[0-9]+( +[0-9]+)* *$/) { nc = NF; for (k = 1; k <= NF; k++) jc[k] = $k; next }   # column header
    if ($0 ~ /^ +[0-9]+ +-?[0-9.]+D[+-][0-9]+/) {
        for (k = 2; k <= NF; k++) { J[jt, $1, jc[k-1]] = dnum($k); if ($1 > jmax) jmax = $1 }
        next }
    inj = 0
}

# ---------------- hyperfine ----------------
/Isotropic Fermi Contact Couplings/ { infc = 1; nfc = 0; getline; next }
infc { if (NF < 6 || $0 ~ /---/) { infc = 0; next }
       nfc++; fc_at[nfc] = $1; fc_iso[nfc] = $2; fc_au[nfc] = $3; fc_mhz[nfc] = $4; fc_g[nfc] = $5; next }

END {
    status = (error ? "error" : (normal ? "normal" : "incomplete"))
    if (mode == "summary") {
        nimag = 0; low = ""
        for (i = 1; i <= nf; i++) { if (fr[i] < 0) nimag++; if (low == "" || fr[i] < low) low = fr[i] }
        e = escf
        h_ev = (homo != "" ? sprintf("%.4f", homo * H2EV) : ""); l_ev = (lumo != "" ? sprintf("%.4f", lumo * H2EV) : "")
        gap = (homo != "" && lumo != "" ? sprintf("%.4f", (lumo - homo) * H2EV) : "")
        if (header) print "name,status,charge,mult,natoms,E_scf,E_post,post_method,ZPE,E_thermal,H,G,Cv,S,HOMO_eV,LUMO_eV,gap_eV,dipole_D,n_imag,lowest_freq,T_K"
        print name "," status "," charge "," mult "," nat "," e "," epost "," lpost "," zpe "," eth "," hh "," gg "," cv "," ss "," h_ev "," l_ev "," gap "," dip "," (nf ? nimag : "") "," low "," T
    }
    else if (mode == "thermo") {
        if (header) print "name,T_K,P_atm,E_scf,ZPE_corr,H_corr,G_corr,E0,H,G,Cv_cal_molK,S_cal_molK"
        for (i = 1; i <= nT; i++)
            print name "," tT[i] "," tP[i] "," tE[i] "," tZ[i] "," tH[i] "," tG[i] "," tE0[i] "," tHH[i] "," tGG[i] "," tCv[i] "," tS[i]
    }
    else if (mode == "charges") {
        if (header) print "name,atom,element,mulliken,apt,npa,hirshfeld,cm5,spin"
        n = nm; if (nn > n) n = nn; if (nh > n) n = nh; if (na > n) n = na
        for (i = 1; i <= n; i++)
            print name "," i "," el[i] "," mul[i] "," apt[i] "," npa[i] "," hir[i] "," cm5[i] "," (spin[i] != "" ? spin[i] : hspin[i])
    }
    else if (mode == "uv") {
        if (header) print "name,state,multiplicity,E_eV,lambda_nm,f,main_transitions(>=10%)"
        for (i = 1; i <= ns; i++) print name "," st[i] "," smult[i] "," sev[i] "," snm[i] "," sf[i] ",\"" strans[i] "\""
    }
    else if (mode == "freq") {
        if (header) print "name,mode,freq_cm1,IR_km_mol,Raman_A4_amu"
        for (i = 1; i <= nf; i++) print name "," i "," fr[i] "," ir[i] "," ra[i]
    }
    else if (mode == "nmr") {
        if (header) print "name,atom,element,shielding_iso_ppm,anisotropy_ppm"
        for (i = 1; i <= nmr_n; i++) print name "," nmr_at[i] "," nmr_el[i] "," nmr_iso[i] "," nmr_an[i]
    }
    else if (mode == "polar") {
        # beta_i = beta_iii + beta_ijj + beta_ikk (Kleinman symmetry), beta_tot = |beta|
        # This is the convention used in most NLO papers. Gaussian's own "x,y,z"
        # vector components are 3 x beta_i, and its "||" value is printed too.
        if (a_iso == "") a_iso = a_iso_simple
        if (header) print "name,alpha_iso_au,alpha_aniso_au,beta_x_au,beta_y_au,beta_z_au,beta_tot_au,beta_tot_1e-30esu,beta_parallel_au(Gaussian)"
        if (!nB) { print name "," a_iso "," a_aniso ",,,,,,"; exit }
        bx = B["xxx"] + B["yxy"] + B["zxz"]; by = B["xxy"] + B["yyy"] + B["zyz"]; bz = B["xxz"] + B["yyz"] + B["zzz"]
        bt = sqrt(bx * bx + by * by + bz * bz)
        printf "%s,%s,%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n", name, a_iso, a_aniso, bx, by, bz, bt, bt * BAU2ESU, bpar
    }
    else if (mode == "hyperfine") {
        if (header) print "name,atom,isotope,fermi_au,fermi_MHz,fermi_Gauss"
        for (i = 1; i <= nfc; i++) print name "," fc_at[i] "," fc_iso[i] "," fc_au[i] "," fc_mhz[i] "," fc_g[i]
    }
    else if (mode == "jcoupling") {
        # pairs="1-2,3-7" limits the output, otherwise all pairs with |J| >= jmin Hz (default 0.5)
        split("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr", S, " ")
        if (jmin == "") jmin = 0.5
        np = split(pairs, pp, ","); for (k = 1; k <= np; k++) { split(pp[k], ab, "-"); want[ab[1] + 0, ab[2] + 0] = 1; want[ab[2] + 0, ab[1] + 0] = 1 }
        if (header) print "name,atom_i,el_i,atom_j,el_j,FC_Hz,SD_Hz,PSO_Hz,DSO_Hz,J_total_Hz"
        for (i = 2; i <= jmax; i++) for (j = 1; j < i; j++) {
            if (!(("J", i, j) in J)) continue
            if (np) { if (!((i, j) in want)) continue } else if (J["J", i, j] < jmin && J["J", i, j] > -jmin) continue
            printf "%s,%d,%s,%d,%s,%.3f,%.3f,%.3f,%.3f,%.3f\n", name, j, S[gz[j]], i, S[gz[i]], J["FC", i, j], J["SD", i, j], J["PSO", i, j], J["DSO", i, j], J["J", i, j]
        }
    }
    else if (mode == "geometry") {
        split("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi", S, " ")
        print ngeo; print name " final geometry (" status ")"
        for (i = 1; i <= ngeo; i++) printf "%-3s %14.8f %14.8f %14.8f\n", S[gz[i]], gx[i], gy[i], gzz[i]
    }
    else { print "unknown mode: " mode > "/dev/stderr"; exit 1 }
}
