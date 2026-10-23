#!/usr/bin/env bash
# =============================================================================
#  redox-table.sh - build the energies.csv used by redox-pka.sh from the summary
#                   tables of any extract script (Gaussian, ORCA, NWChem)
#
#  usage: redox-table.sh -g gas.csv -s solvent.csv [-H high_gas.csv] [-z 0.9806] > energies.csv
#    -g  summary of the gas phase opt+freq jobs (gives G and E at the low level)
#    -s  summary of the single points in implicit solvent (E_solv)
#    -H  summary of gas phase single points at the higher level used for -s;
#        then G_gas = E_high + (G_low - E_low), the usual composite recipe, and
#        dG_solv = E_solv - E_high. Without -H, E_gas is the low level energy.
#    -z  scale the thermal free energy correction (G - E) by this factor
#        (the original workflow used 0.9806 for B3LYP/6-31G(d) zero point energies)
#  Rows are matched by molecule name; the row that has G is taken from -g.
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"
gas=""; sol=""; high=""; z=1
while getopts "g:s:H:z:h" o; do case $o in g) gas=$OPTARG ;; s) sol=$OPTARG ;; H) high=$OPTARG ;; z) z=$OPTARG ;; *) sed -n '2,15p' "$0"; exit 0 ;; esac; done
[[ -f "$gas" && -f "$sol" ]] || die "need -g gas.csv and -s solvent.csv"
pick() {   # pick FILE -> name E [G]
    awk -F, 'NR==1 { for (i=1;i<=NF;i++) c[$i]=i; e = c["E_scf"] ? c["E_scf"] : (c["E_final"] ? c["E_final"] : c["E_total"]); next }
             { g = c["G"] ? $c["G"] : ""; if (!($1 in s) || g != "") { s[$1] = 1; E[$1] = $e; G[$1] = g } }
             END { for (k in E) print k, E[k], G[k] }' "$1"
}
echo "name,G_gas,E_gas,E_solv"
awk -v z="$z" 'FILENAME == ARGV[1] { Eg[$1] = $2; Gg[$1] = $3; next }
               FILENAME == ARGV[2] { Es[$1] = $2; next }
               { Eh[$1] = $2 }
               END { for (k in Eg) { if (Gg[k] == "" || !(k in Es)) { print "warning: " k " lacks G or a solvent energy" > "/dev/stderr"; continue }
                       e = (k in Eh) ? Eh[k] : Eg[k]
                       printf "%s,%.8f,%.8f,%.8f\n", k, e + z * (Gg[k] - Eg[k]), e, Es[k] } }' \
    <(pick "$gas") <(pick "$sol") <( [[ -n "$high" ]] && pick "$high" ) | sort
