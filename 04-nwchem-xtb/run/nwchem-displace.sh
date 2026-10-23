#!/usr/bin/env bash
# =============================================================================
#  nwchem-displace.sh - follow the imaginary mode of a finished NWChem frequency job
#
#  NWChem prints the structure displaced by +/-100 % along each imaginary mode.
#  * a "minimum" with an imaginary frequency: re-optimise from NAME-minus.xyz
#    (or NAME-plus.xyz); it lies downhill along the offending mode
#  * a transition state: optimise both NAME-plus and NAME-minus; they should
#    fall into the reactant and the product, a quick check that the TS
#    connects the intended minima (a poor man's IRC)
#
#  usage: nwchem-displace.sh [-m MODE] [-p PERCENT] LOG [OUTDIR]
#    -m  mode number (default 1, the most negative)   -p 100 or 50 (default 100)
#  writes OUTDIR/NAME-plus.xyz and OUTDIR/NAME-minus.xyz (OUTDIR default: restart)
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"
mode=1; pct=100
while getopts "m:p:h" o; do case $o in m) mode=$OPTARG ;; p) pct=$OPTARG ;; *) sed -n '2,16p' "$0"; exit 0 ;; esac; done
shift $((OPTIND-1))
log=${1:?give an NWChem log}; out=${2:-restart}; mkdir -p "$out"
name=$(strip_ext "$log"); name=${name%-opt}; name=${name%-ts}
for sign in plus minus; do
    p=$pct; [[ $sign == minus ]] && p=-$pct
    awk -v want="$(printf '%.1f%%' "$p")" -v m="$mode" -v nm="$name-$sign" '
        $1 == "Geometry" && $2 == "after" && $3 == want && $7 == m";" { n = 0; getline; grab = 1; next }
        grab { if (NF < 6) { grab = 0; next } n++; e[n] = $2; x[n] = $4; y[n] = $5; z[n] = $6 }
        END { if (!n) exit 1; print n; print nm " displaced along mode " m
              for (i = 1; i <= n; i++) printf "%-3s %14.8f %14.8f %14.8f\n", e[i], x[i], y[i], z[i] }' "$log" > "$out/$name-$sign.xyz" ||
        { rm -f "$out/$name-$sign.xyz"; die "no displaced geometry for mode $mode in $log (no imaginary frequency?)"; }
done
info "$out/$name-plus.xyz and $out/$name-minus.xyz written"
