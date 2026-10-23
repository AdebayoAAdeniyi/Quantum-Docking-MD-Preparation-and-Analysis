#!/usr/bin/env bash
# =============================================================================
#  gamess-extract.sh - GAMESS-US / Firefly results of a batch as CSV tables
#
#  usage: gamess-extract.sh -k KIND [-l list] [-p parent] [-o outdir]
#    KIND: summary | thermo | freq | charges | neda | geometry | all
#  Logs are found as PARENT/NAME/NAME-*.log (default parent gamess-jobs).
#  "neda" also writes neda-wide.csv: one row per molecule, one column per
#  energy component, ready for plotting or correlation with experiment.
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"
AWKM="$QDM_ROOT/lib/awk/gamess.awk"
kind=""; list=""; parent=gamess-jobs; out=gamess-results
while getopts "k:l:p:o:h" o; do
    case $o in k) kind=$OPTARG ;; l) list=$OPTARG ;; p) parent=$OPTARG ;; o) out=$OPTARG ;; *) sed -n '2,10p' "$0"; exit 0 ;; esac
done
[[ -n "$kind" ]] || die "choose -k (see -h)"
if [[ -n "$list" ]]; then mapfile -t names < <(read_list "$list")
else mapfile -t names < <(find "$parent" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort); fi
logs=(); for n in "${names[@]}"; do for f in "$parent/$n/$n"-*.log; do [[ -f "$f" ]] && logs+=("$f"); done; done
((${#logs[@]})) || die "no GAMESS/Firefly logs found in $parent"
mkdir -p "$out"
table() {
    local k=$1 csv="$out/$1.csv" first=1 f n s
    : > "$csv"
    for f in "${logs[@]}"; do
        n=$(basename "$(dirname "$f")"); s=$(strip_ext "$f"); s=${s#"$n"-}
        awk -v mode="$k" -v name="$n,$s" -v header=$first -f "$AWKM" "$f" |
            { if ((first)); then sed '1s/^name,/name,step,/'; else cat; fi; } >> "$csv"
        [[ -s "$csv" ]] && first=0
    done
    if [[ $k != summary && $(wc -l < "$csv") -le 1 ]]; then warn "no $k data"; rm -f "$csv"; return; fi
    info "$k: $csv"
}
neda_wide() {
    [[ -f "$out/neda.csv" ]] || return 0
    awk -F, 'NR>1 { k=$1; if (!(k in seen)) { seen[k]=1; order[++n]=k } if (!($3 in c)) { c[$3]=1; col[++m]=$3 } v[k,$3]=$4 }
        END { printf "name"; for (j=1;j<=m;j++) printf ",%s", col[j]; print ""
              for (i=1;i<=n;i++) { printf "%s", order[i]; for (j=1;j<=m;j++) printf ",%s", v[order[i],col[j]]; print "" } }' \
        "$out/neda.csv" > "$out/neda-wide.csv"
    info "neda: $out/neda-wide.csv"
}
geometry() {
    for f in "${logs[@]}"; do
        grep -q "COORDINATES OF ALL ATOMS ARE (ANGS)" "$f" || continue
        awk -v mode=geometry -v name="$(strip_ext "$f")" -f "$AWKM" "$f" > "$out/$(strip_ext "$f")-final.xyz"
    done
}
case $kind in
    all) for k in summary thermo freq charges neda; do table $k; done; neda_wide; geometry ;;
    neda) table neda; neda_wide ;;
    geometry) geometry ;;
    summary|thermo|freq|charges) table "$kind" ;;
    *) die "unknown kind $kind" ;;
esac
