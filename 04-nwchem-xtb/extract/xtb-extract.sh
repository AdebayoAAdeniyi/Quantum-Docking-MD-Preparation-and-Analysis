#!/usr/bin/env bash
# =============================================================================
#  xtb-extract.sh - one table for a whole xtb screening run
#  usage: xtb-extract.sh [-l list] [-p parent] [-o outdir]
#  writes xtb-screen.csv: energy, H, G, gap, dipole, imaginary modes, IP, EA and
#  the global descriptors mu = -(IP+EA)/2, eta = (IP-EA)/2, omega = mu^2/(2 eta),
#  plus xtb-fukui.csv with the condensed Fukui indices of every atom
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"
A="$QDM_ROOT/lib/awk/xtb.awk"; list=""; parent=xtb-jobs; out=xtb-results
while getopts "l:p:o:h" o; do case $o in l) list=$OPTARG ;; p) parent=$OPTARG ;; o) out=$OPTARG ;; *) sed -n '2,8p' "$0"; exit 0 ;; esac; done
if [[ -n "$list" ]]; then mapfile -t names < <(read_list "$list"); else mapfile -t names < <(find "$parent" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort); fi
mkdir -p "$out"; csv="$out/xtb-screen.csv"; fk="$out/xtb-fukui.csv"
echo "name,status,E_total,H,G,gap_eV,dipole_D,n_imag,IP_eV,EA_eV,mu_eV,eta_eV,omega_eV" > "$csv"; : > "$fk"; first=1
for n in "${names[@]}"; do
    d="$parent/$n"; o=$(ls "$d/$n-ohess.out" "$d/$n-opt.out" 2>/dev/null | head -1)
    [[ -n "$o" ]] || { warn "$n: no xtb optimisation output"; continue; }
    IFS=, read -r _ st _ e h g _ _ _ gap dip nim <<< "$(awk -v mode=summary -f "$A" "$o")"
    ip=""; ea=""; [[ -f "$d/$n-vipea.out" ]] && IFS=, read -r _ ip ea <<< "$(awk -v mode=vipea -f "$A" "$d/$n-vipea.out")"
    awk -v n="$n" -v st="$st" -v e="$e" -v h="$h" -v g="$g" -v gap="$gap" -v dip="$dip" -v nim="$nim" -v ip="$ip" -v ea="$ea" 'BEGIN{
        d = ""; if (ip != "" && ea != "") { mu = -(ip + ea) / 2; eta = (ip - ea) / 2; d = sprintf("%.4f,%.4f,%.4f", mu, eta, mu * mu / (2 * eta)) } else d = ",,"
        printf "%s,%s,%s,%s,%s,%.4f,%s,%s,%s,%s,%s\n", n, st, e, h, g, gap, dip, nim, ip, ea, d }' >> "$csv"
    [[ -f "$d/$n-vfukui.out" ]] && { awk -v mode=fukui -v name="$n" -v header=$first -f "$A" "$d/$n-vfukui.out" >> "$fk"; first=0; }
done
[[ -s "$fk" ]] || rm -f "$fk"
info "xtb: $csv $( [[ -f "$fk" ]] && echo "and $fk")"
