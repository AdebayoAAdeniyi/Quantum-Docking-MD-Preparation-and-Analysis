#!/usr/bin/env bash
# =============================================================================
#  nwchem-extract.sh - NWChem results of a batch as CSV tables
#  usage: nwchem-extract.sh -k KIND [-l list] [-p parent] [-o outdir]
#    KIND: summary | thermo | freq | neb | bsse | geometry | all
#  Logs: PARENT/NAME/NAME-*.log (default parent nwchem-jobs)
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"; source "$QDM_ROOT/lib/table.sh"
A="$QDM_ROOT/lib/awk/nwchem.awk"; kind=""; list=""; parent=nwchem-jobs; out=nwchem-results
while getopts "k:l:p:o:h" o; do case $o in k) kind=$OPTARG ;; l) list=$OPTARG ;; p) parent=$OPTARG ;; o) out=$OPTARG ;; *) sed -n '2,6p' "$0"; exit 0 ;; esac; done
[[ -n "$kind" ]] || die "choose -k"
if [[ -n "$list" ]]; then mapfile -t names < <(read_list "$list"); else mapfile -t names < <(find "$parent" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort); fi
logs=(); for n in "${names[@]}"; do for f in "$parent/$n/$n"-*.log; do [[ -f "$f" ]] && logs+=("$f"); done; done
((${#logs[@]})) || die "no NWChem logs in $parent"; mkdir -p "$out"
geometry() { for f in "${logs[@]}"; do grep -q "Output coordinates in angstroms" "$f" && awk -v mode=geometry -v name="$(strip_ext "$f")" -f "$A" "$f" > "$out/$(strip_ext "$f")-final.xyz"; done; info "geometry: $out/*-final.xyz"; }
case $kind in
    all) for k in summary thermo freq neb bsse; do
             sel=("${logs[@]}"); [[ $k == neb ]] && sel=($(printf '%s\n' "${logs[@]}" | grep -- '-neb\.log$'))
             [[ $k == bsse ]] && sel=($(printf '%s\n' "${logs[@]}" | grep -- '-bsse\.log$'))
             ((${#sel[@]})) && build_table "$A" $k "$out/$k.csv" "${sel[@]}"; done; geometry ;;
    geometry) geometry ;;
    summary|thermo|freq|neb|bsse) build_table "$A" "$kind" "$out/$kind.csv" "${logs[@]}" ;;
    *) die "unknown kind $kind" ;;
esac
