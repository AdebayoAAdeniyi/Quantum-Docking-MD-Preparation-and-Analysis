# =============================================================================
#  xtb.awk - read xtb 6.x output
#  modes: summary  name,status,charge,E_total,H,G,ZPE,HOMO_eV,LUMO_eV,gap_eV,dipole_D,n_imag
#         vipea    name,IP_eV,EA_eV          (vertical, delta SCC)
#         fukui    name,atom,element,f+,f-,f0  (as printed by xtb)
#         geometry not available here (use NAME-xtbopt.xyz)
# =============================================================================
BEGIN { if (name == "") name = "mol" }
/normal termination of xtb/ { normal = 1 }
/\[ERROR\]|abnormal termination/ { err = 1 }
/:: total charge +/ { chg = sprintf("%.0f", $4) + 0 }
/\(HOMO\)$/ { homo = $(NF - 1) }
/\(LUMO\)$/ { lumo = $(NF - 1) }
/\| TOTAL ENERGY / { e = $4 }
/\| TOTAL ENTHALPY / { h = $4 }
/\| TOTAL FREE ENERGY / { g = $5 }
/\| HOMO-LUMO GAP / { gap = $4 }
/:: zero point energy / { zpe = $5 }
/^molecular dipole:/ { indip = 1 }
indip && /^ +full:/ { dip = $NF; indip = 0 }
/# imaginary freq\./ { nim = $(NF - 1) }
/delta SCC IP \(eV\):/ { ip = $NF }
/delta SCC EA \(eV\):/ { ea = $NF }
/ +# +f\(\+\) +f\(-\) +f\(0\)/ { inf = 1; nf = 0; next }
inf { if ($0 !~ /^ +[0-9]+[A-Za-z]+ +-?[0-9.]+/) { inf = 0; next }
      nf++; at = $1; el = at; gsub(/[0-9]/, "", el); gsub(/[A-Za-z]/, "", at); fa[nf] = at; fe[nf] = el; fp[nf] = $2; fm[nf] = $3; f0[nf] = $4 }
END {
    st = (err ? "error" : (normal ? "normal" : "incomplete"))
    if (mode == "summary") {
        if (header) print "name,status,charge,E_total,H,G,ZPE,HOMO_eV,LUMO_eV,gap_eV,dipole_D,n_imag"
        print name "," st "," chg "," e "," h "," g "," zpe "," homo "," lumo "," gap "," dip "," nim
    } else if (mode == "vipea") {
        if (header) print "name,IP_eV,EA_eV"
        print name "," ip "," ea
    } else if (mode == "fukui") {
        if (header) print "name,atom,element,f_plus,f_minus,f_zero"
        for (i = 1; i <= nf; i++) print name "," fa[i] "," fe[i] "," fp[i] "," fm[i] "," f0[i]
    } else { print "unknown mode " mode > "/dev/stderr"; exit 1 }
}
