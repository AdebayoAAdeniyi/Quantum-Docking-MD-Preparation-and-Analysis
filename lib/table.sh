#!/usr/bin/env bash
# build_table AWKFILE KIND CSV FILE... : run one parser mode over many outputs and
# write one CSV with "name,step" in front (name = job folder, step = file suffix)
build_table() {
    local awkf=$1 k=$2 csv=$3 first=1 f n s; shift 3
    : > "$csv"
    for f in "$@"; do
        n=$(basename "$(dirname "$f")"); s=$(strip_ext "$f"); s=${s#"$n"}; s=${s#-}; s=${s:-main}
        awk -v mode="$k" -v name="$n,$s" -v header=$first -f "$awkf" "$f" 2>/dev/null |
            { if ((first)); then sed '1s/^name,/name,step,/'; else cat; fi; } >> "$csv"
        [[ -s "$csv" ]] && first=0
    done
    if [[ $k != summary && $(wc -l < "$csv") -le 1 ]]; then warn "no $k data"; rm -f "$csv"; return 0; fi
    info "$k: $csv ($(($(wc -l < "$csv") - 1)) rows)"
}
