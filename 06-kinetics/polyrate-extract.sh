#!/usr/bin/env bash
# =============================================================================
#  polyrate-extract.sh - rate constants from POLYRATE long output (.fu6) files
#
#  usage: polyrate-extract.sh [-l list] [-p parent] [-o outdir]
#    reads PARENT/NAME/NAME.fu6 and copies the final table that starts with
#    "T(K)  Forward  Reverse" (the rate constants of the last method printed,
#    e.g. CVT/SCT) to polyrate-rates.csv, then fits ln k = ln A - Ea/RT.
#  Same table that the original extract-polyrate script read; no POLYRATE output
#  was available to test it here.
# =============================================================================
source "$(dirname "$0")/../lib/common.sh"
list=""; parent=polyrate-jobs; out=polyrate-results
while getopts "l:p:o:h" o; do case $o in l) list=$OPTARG ;; p) parent=$OPTARG ;; o) out=$OPTARG ;; *) sed -n '2,11p' "$0"; exit 0 ;; esac; done
if [[ -n "$list" ]]; then mapfile -t names < <(read_list "$list"); else mapfile -t names < <(find "$parent" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort); fi
mkdir -p "$out"; csv="$out/polyrate-rates.csv"; echo "name,T_K,k_forward,k_reverse" > "$csv"
for n in "${names[@]}"; do
    f="$parent/$n/$n.fu6"; [[ -f "$f" ]] || { warn "$n: no $f"; continue; }
    awk -v n="$n" '{ sub(/\r$/, "") }
        /T\(K\) +Forward +Reverse/ { nr = 0; delete T; intab = 1; next }
        intab && /^ *[0-9.]+ +[0-9.]+[EeDd][-+][0-9]+/ { nr++; T[nr] = $1; F[nr] = $2; R[nr] = $3; next }
        intab && nr && !/^ *$/ { intab = 0 }
        END { for (i = 1; i <= nr; i++) { gsub(/[dD]/, "E", F[i]); gsub(/[dD]/, "E", R[i]); print n "," T[i] "," F[i] "," R[i] } }' "$f" >> "$csv"
done
awk -F, 'NR>1 && $3 > 0 { x = 1/$2; y = log($3); k = $1; n[k]++; sx[k]+=x; sy[k]+=y; sxx[k]+=x*x; sxy[k]+=x*y }
    END { print "name,n_T,Ea_kJ_mol,A"; for (k in n) if (n[k] > 2) { b = (n[k]*sxy[k]-sx[k]*sy[k])/(n[k]*sxx[k]-sx[k]^2); a = (sy[k]-b*sx[k])/n[k]
          printf "%s,%d,%.2f,%.4e\n", k, n[k], -b*8.314462618/1000, exp(a) } }' "$csv" > "$out/polyrate-arrhenius.csv"
info "POLYRATE: $csv, $out/polyrate-arrhenius.csv"
