#!/usr/bin/env bash
# =============================================================================
#  plot-profile.sh - free energy profiles (reactants -> TS -> products) for every
#                    reaction in kinetics-results/profile.dat, drawn as the usual
#                    horizontal levels joined by dashed lines
#  usage: plot-profile.sh [-i kinetics-results/profile.dat] [-o profile.png]
# =============================================================================
source "$(dirname "$0")/../lib/common.sh"; require gnuplot
in=kinetics-results/profile.dat; out=profile.png
while getopts "i:o:h" o; do case $o in i) in=$OPTARG ;; o) out=$OPTARG ;; *) sed -n '2,7p' "$0"; exit 0 ;; esac; done
[[ -s "$in" ]] || die "no $in (run kinetics.sh first)"
base=${out%.*}-levels
# each level becomes a short horizontal segment; consecutive levels are joined
awk '{ if ($1 != last && NR > 1) print "\n"; last = $1; printf "%s %.3f %s\n%s %.3f %s\n", $1, $2*3, $3, $1, $2*3+1.4, $3 }' "$in" > "$base.dat"
names=$(awk '!s[$1]++{print $1}' "$in")
plots=(); i=0
for n in $names; do
    plots+=("'$base.dat' index $i u 2:3 w l lw 4 t '$n'" "'$base.dat' index $i every 2::1 u (\$2-0.7):3:(sprintf('%.1f',\$3)) w labels offset 0,1 notitle")
    i=$((i + 1))
done
cat > "$base.gp" <<EOG
set terminal pngcairo size 900,650 enhanced font 'Arial,13'
set output '$out'
set ylabel '{/Symbol D}G (kcal mol^{-1})'
set xtics ('reactants' 0.7, 'TS' 3.7, 'products' 6.7)
set xrange [-0.5:7.9]
set grid ytics lc rgb '#dddddd'
set key top right
plot $(IFS=,; echo "${plots[*]}")
EOG
gnuplot "$base.gp" && info "profile: $out"
