#!/usr/bin/env bash
# =============================================================================
#  rosetta-scores.sh - Rosetta score files to CSV, ranked
#
#  Reads .sc, .fasc and silent score tables, writes one CSV with every score
#  term and a second sorted on the term you choose.
#
#  usage: rosetta-scores.sh [options] [FILE_OR_DIR ...]
#    -o DIR     output directory (default: scores)
#    -k TERM    term to sort on (default: total_score, use I_sc for interfaces)
#    -r         sort descending (default is ascending, best first for energies)
#    -n N       print the best N rows (default: 5)
#    -h         this help
#
#  Why this is not a sort command
#  ------------------------------
#  A Rosetta score file carries a header row, and the original scripts ranked
#  models with "cat *.fasc | sort -nk2" and "sort -nk5", then took the first
#  line. Three things go wrong with that:
#
#  1. Concatenating several score files repeats the header, and a header line
#     sorts numerically as zero. With favourable (negative) scores it lands in
#     the middle; with positive scores it lands first, and the "best model"
#     read off the top is the word "description".
#  2. Column 2 is only total_score, and column 5 only I_sc, for one particular
#     set of score terms. Change the weights or the protocol and the columns
#     move, silently ranking on whatever now sits there.
#  3. Score files from different protocols have different columns, so pooling
#     them puts unlike numbers in one sort.
#
#  Terms are looked up by name from each file's own header here, header rows
#  are dropped, and files whose columns disagree are reported rather than
#  merged. Both the whitespace "SCORE:" layout and the newer one-JSON-object
#  per line .fasc layout are read.
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../../lib/common.sh"

out="scores"; key="total_score"; desc=0; top=5
while getopts "o:k:rn:h" o; do case $o in
    o) out=$OPTARG ;; k) key=$OPTARG ;; r) desc=1 ;; n) top=$OPTARG ;;
    h) sed -n '2,/^# ====/p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) die "unknown option; -h for usage" ;;
esac; done
shift $((OPTIND - 1))
(($#)) || set -- .

shopt -s nullglob
files=()
for a in "$@"; do
    if [[ -d "$a" ]]; then files+=("$a"/*.fasc "$a"/*.sc)
    elif [[ -f "$a" ]]; then files+=("$a")
    else warn "not found: $a"; fi
done
(( ${#files[@]} )) || die "no .sc or .fasc score files found"

mkdir -p "$out"
info "reading ${#files[@]} score file(s)"

# ---------------------------------------------------------------------------
# Parse every file into a long table: file, description, term, value
# ---------------------------------------------------------------------------
awk '
    function trim(s) { gsub(/^[[:space:]"]+|[[:space:]",]+$/, "", s); return s }

    FNR == 1 { nh = 0; delete H }

    # ---- newer .fasc layout: one JSON object per line ----
    /^[[:space:]]*\{/ {
        line = $0
        gsub(/^[[:space:]]*\{|\}[[:space:]]*$/, "", line)
        n = split(line, kv, ",")
        d = ""; delete V; nv = 0
        for (i = 1; i <= n; i++) {
            p = index(kv[i], ":")
            if (!p) continue
            k = trim(substr(kv[i], 1, p - 1))
            v = trim(substr(kv[i], p + 1))
            if (k == "decoy" || k == "description") d = v
            else { nv++; K[nv] = k; Vv[nv] = v }
        }
        if (d == "") d = sprintf("%s:%d", FILENAME, FNR)
        for (i = 1; i <= nv; i++) printf "%s\t%s\t%s\t%s\n", FILENAME, d, K[i], Vv[i]
        next
    }

    # ---- classic layout: a SCORE: header row then SCORE: data rows ----
    $1 == "SCORE:" || $1 == "SCORE" {
        # A header row is the one whose second field is not a number.
        if ($2 + 0 == 0 && $2 != "0" && $2 !~ /^-?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/) {
            nh = NF; for (i = 2; i <= NF; i++) H[i] = $i
            next
        }
        if (nh == 0) { bad[FILENAME] = 1; next }    # data before any header
        d = $NF
        for (i = 2; i <= NF && i <= nh; i++) {
            if (H[i] == "description") continue
            printf "%s\t%s\t%s\t%s\n", FILENAME, d, H[i], $i
        }
        next
    }

    END { for (f in bad) printf "%s has score rows but no header row\n", f > "/dev/stderr" }
' "${files[@]}" > "$out/long.tsv"

grep -q . "$out/long.tsv" || die "no score rows parsed; are these Rosetta score files?"

# ---------------------------------------------------------------------------
# Long table to a wide CSV, with the union of terms as columns
# ---------------------------------------------------------------------------
awk -F'\t' '
    { V[$2 "\t" $3] = $4; if (!($2 in seenD)) { seenD[$2] = 1; D[++nd] = $2; SRC[$2] = $1 }
      if (!($3 in seenT)) { seenT[$3] = 1; T[++nt] = $3 } }
    END {
        printf "description,source"
        for (j = 1; j <= nt; j++) printf ",%s", T[j]
        printf "\n"
        for (i = 1; i <= nd; i++) {
            printf "%s,%s", D[i], SRC[D[i]]
            for (j = 1; j <= nt; j++) {
                k = D[i] "\t" T[j]
                printf ",%s", (k in V) ? V[k] : ""
            }
            printf "\n"
        }
    }
' "$out/long.tsv" > "$out/scores.csv"

n_models=$(( $(wc -l < "$out/scores.csv") - 1 ))
n_terms=$(head -1 "$out/scores.csv" | awk -F, '{print NF - 2}')
info "wrote $out/scores.csv: $n_models model(s), $n_terms score term(s)"

# Report terms that are missing for some models, which is the sign that score
# files from different protocols have been pooled.
awk -F, 'NR==1{for(i=3;i<=NF;i++)h[i]=$i; next}
         {for(i=3;i<=NF;i++) if($i=="") miss[h[i]]++}
         END{for(t in miss) printf "  %s missing for %d model(s)\n", t, miss[t]}' \
    "$out/scores.csv" | sort > "$out/incomplete-terms.txt"
if [[ -s "$out/incomplete-terms.txt" ]]; then
    warn "some score terms are not present for every model:"
    cat "$out/incomplete-terms.txt" >&2
    warn "that usually means score files from different protocols were pooled"
fi

# ---------------------------------------------------------------------------
# Rank on the chosen term, by name
# ---------------------------------------------------------------------------
col=$(head -1 "$out/scores.csv" | awk -F, -v k="$key" '{for(i=1;i<=NF;i++) if($i==k){print i; exit}}')
if [[ -z "$col" ]]; then
    warn "term '$key' is not in these score files. Available:"
    head -1 "$out/scores.csv" | tr ',' '\n' | sed 's/^/  /' | tail -n +3 >&2
    die "choose one with -k"
fi

{
    head -1 "$out/scores.csv"
    tail -n +2 "$out/scores.csv" | awk -F, -v c="$col" '$c != ""' |
        sort -t, -g $( ((desc)) && echo -r ) -k"$col,$col"
} > "$out/scores-sorted.csv"

blank=$(tail -n +2 "$out/scores.csv" | awk -F, -v c="$col" '$c == ""' | wc -l)
(( blank )) && warn "$blank model(s) have no value for '$key' and were left out of the ranking"

info "wrote $out/scores-sorted.csv (ranked on $key, $( ((desc)) && echo descending || echo ascending ))"
printf '\nBest %s by %s:\n' "$top" "$key"
{ head -1 "$out/scores-sorted.csv"; tail -n +2 "$out/scores-sorted.csv" | head -"$top"; } |
    awk -F, -v c="$col" 'NR==1{printf "  %-45s %14s\n","model",$c; next}
                         {printf "  %-45s %14s\n",$1,$c}'
