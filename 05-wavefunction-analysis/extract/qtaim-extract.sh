#!/usr/bin/env bash
# =============================================================================
#  qtaim-extract.sh - bond critical point table for a batch of AIMAll results
#  usage: qtaim-extract.sh [-l list] [-p parent] [-o outdir] [-e "Fe,Ru"] [-r 0.04]
#    writes bcp-all.csv (every BCP with H, -G/V and bond class),
#    bcp-weak.csv  (rho <= -r a.u., the hydrogen-bond-like and weak contacts; default 0.04),
#    bcp-metal.csv (BCPs that involve the elements given by -e)
#  The two filters reproduce the selections of the original select-bonds-qtaim.
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"
list=""; parent=gaussian-jobs; out=qtaim-results; els=""; rmax=0.04
while getopts "l:p:o:e:r:h" o; do case $o in l) list=$OPTARG ;; p) parent=$OPTARG ;; o) out=$OPTARG ;; e) els=$OPTARG ;; r) rmax=$OPTARG ;; *) sed -n '2,10p' "$0"; exit 0 ;; esac; done
if [[ -n "$list" ]]; then mapfile -t names < <(read_list "$list"); else mapfile -t names < <(find "$parent" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort); fi
mkdir -p "$out"; all="$out/bcp-all.csv"; first=1; : > "$all"
for n in "${names[@]}"; do
    for f in "$parent/$n"/*.sum "$parent/$n"/*.mgpviz; do
        [[ -f "$f" ]] || continue
        awk -v name="$(strip_ext "$f")" -v header=$first -f "$QDM_ROOT/lib/awk/aimall.awk" "$f" >> "$all"; first=0
    done
done
[[ -s "$all" ]] || die "no .sum or .mgpviz files found under $parent"
awk -F, -v r="$rmax" 'NR==1 || $5 <= r' "$all" > "$out/bcp-weak.csv"
if [[ -n "$els" ]]; then
    awk -F, -v e="$els" 'BEGIN{n=split(e,E,",")} NR==1{print; next} { a=$3; b=$4; gsub(/[0-9]+/,"",a); gsub(/[0-9]+/,"",b); for(i=1;i<=n;i++) if (a==E[i] || b==E[i]) { print; break } }' "$all" > "$out/bcp-metal.csv"
fi
info "QTAIM: $all ($(($(wc -l < "$all") - 1)) BCPs), $out/bcp-weak.csv$( [[ -n "$els" ]] && echo ", $out/bcp-metal.csv")"
