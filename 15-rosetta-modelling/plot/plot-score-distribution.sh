#!/usr/bin/env bash
# =============================================================================
#  plot-score-distribution.sh - score funnel and score histogram
#
#  usage: plot-score-distribution.sh -i scores.csv [options]
#    -i FILE   CSV from extract/rosetta-scores.sh
#    -o DIR    output directory (default: plots)
#    -y TERM   score on the y axis (default: total_score)
#    -x TERM   quantity on the x axis for the funnel (default: rms)
#    -h        this help
#
#  Two plots, because they answer different questions.
#
#  The histogram shows the spread of scores over the models generated. A single
#  narrow peak means the sampling converged on one arrangement; a broad or
#  multi-peaked distribution means it did not, and the lowest-scoring model is
#  then one sample from a wide distribution rather than a converged answer.
#
#  The funnel plots score against a distance measure, usually RMSD to a
#  reference. A real funnel narrows towards low RMSD, which is what says the
#  score function can tell near-native from decoy. A flat cloud says it cannot,
#  and in that case the best-scoring model is not evidence of anything. Without
#  a reference structure there is no RMSD column and only the histogram is drawn.
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../../lib/common.sh"

inp=""; out="plots"; yterm="total_score"; xterm="rms"
while getopts "i:o:y:x:h" o; do case $o in
    i) inp=$OPTARG ;; o) out=$OPTARG ;; y) yterm=$OPTARG ;; x) xterm=$OPTARG ;;
    h) sed -n '2,/^# ====/p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) die "unknown option; -h for usage" ;;
esac; done

[[ -n "$inp" ]] || die "a score CSV is required (-i); -h for usage"
[[ -f "$inp" ]] || die "not found: $inp"
require gnuplot
mkdir -p "$out"

col_of() { head -1 "$1" | awk -F, -v k="$2" '{for(i=1;i<=NF;i++) if($i==k){print i; exit}}'; }
ycol=$(col_of "$inp" "$yterm")
[[ -n "$ycol" ]] || die "term '$yterm' not in $inp; columns: $(head -1 "$inp" | tr ',' ' ')"
xcol=$(col_of "$inp" "$xterm")

tail -n +2 "$inp" | awk -F, -v y="$ycol" '$y != "" {print $y}' > "$out/y.dat"
n=$(wc -l < "$out/y.dat")
(( n > 0 )) || die "no values for '$yterm'"

read -r ymin ymax ymean <<<"$(awk '{s+=$1; if(NR==1||$1<lo)lo=$1; if(NR==1||$1>hi)hi=$1}
                                  END{printf "%.3f %.3f %.3f", lo, hi, s/NR}' "$out/y.dat")"
info "$n model(s), $yterm from $ymin to $ymax, mean $ymean"

gnuplot <<EOF
set terminal pngcairo size 900,600 enhanced font 'Helvetica,12'
set output '$out/score-histogram.png'
set title '$yterm over $n model(s)'
set xlabel '$yterm'
set ylabel 'models'
set style fill solid 0.5
set boxwidth 0.9 relative
binw = ($ymax - $ymin) / 40.0
if (binw <= 0) { binw = 1.0 }
bin(x) = binw * floor(x / binw) + binw / 2.0
set grid ytics
plot '$out/y.dat' using (bin(\$1)):(1.0) smooth freq with boxes lc rgb '#4878a8' notitle
EOF
info "wrote $out/score-histogram.png"

if [[ -n "$xcol" ]]; then
    tail -n +2 "$inp" | awk -F, -v x="$xcol" -v y="$ycol" '$x != "" && $y != "" {print $x, $y}' > "$out/funnel.dat"
    if [[ -s "$out/funnel.dat" ]]; then
        gnuplot <<EOF
set terminal pngcairo size 900,600 enhanced font 'Helvetica,12'
set output '$out/score-funnel.png'
set title 'Score funnel: $yterm against $xterm'
set xlabel '$xterm'
set ylabel '$yterm'
set grid
plot '$out/funnel.dat' using 1:2 with points pt 7 ps 0.5 lc rgb '#a84848' notitle
EOF
        info "wrote $out/score-funnel.png"
    fi
else
    info "no '$xterm' column, so no funnel plot. Run with -in:file:native to get RMSD."
fi
rm -f "$out/y.dat"
