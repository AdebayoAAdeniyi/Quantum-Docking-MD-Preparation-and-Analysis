#!/usr/bin/env bash
# =============================================================================
#  plot-neb.sh - reaction path energy profiles from neb.csv (orca-extract.sh -k neb)
#
#  usage: plot-neb.sh [-i neb.csv] [-n names] [-o neb-profile.png]
#  One curve per reaction, relative energies in kcal/mol, the climbing image
#  (CI) and the refined transition state (TS) are labelled.
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"
require gnuplot
csv=orca-results/neb.csv; names=""; out=neb-profile.png
while getopts "i:n:o:h" o; do
    case $o in i) csv=$OPTARG ;; n) names=$OPTARG ;; o) out=$OPTARG ;; *) sed -n '2,8p' "$0"; exit 0 ;; esac
done
[[ -f "$csv" ]] || die "no table $csv"
[[ -z "$names" ]] && names=$(awk -F, 'NR>1 && !s[$1]++ {print $1}' "$csv" | paste -sd,)
base=${out%.*}; plots=()
for n in ${names//,/ }; do
    # images in path order; the TS row (if any) is drawn as a point at the CI position
    awk -F, -v n="$n" 'NR>1 && $1==n && $3!="TS" {print $3, $5, $7}' "$csv" > "$base-$n.dat"
    awk -F, -v n="$n" 'NR>1 && $1==n && $7=="CI" {ci=$3} NR>1 && $1==n && $3=="TS" {ts=$5} END{if (ts!="") print ci, ts, "TS"}' "$csv" > "$base-$n-ts.dat"
    plots+=("'$base-$n.dat' u 1:2 w lp lw 2 pt 7 t '$n'" "'$base-$n.dat' u 1:2:3 w labels offset 0,1 notitle")
    [[ -s "$base-$n-ts.dat" ]] && plots+=("'$base-$n-ts.dat' u 1:2 w p pt 5 ps 1.5 notitle" "'$base-$n-ts.dat' u 1:2:3 w labels offset 1.5,0 notitle")
done
cat > "$base.gp" <<EOG
set terminal pngcairo size 1000,700 enhanced font 'Arial,13'
set output '$out'
set xlabel 'NEB image'
set ylabel '{/Symbol D}E (kcal mol^{-1})'
set grid lc rgb '#dddddd'
set key top left
plot $(IFS=,; echo "${plots[*]}")
EOG
gnuplot "$base.gp" && info "NEB profile: $out"
