#!/usr/bin/env bash
# =============================================================================
#  lammps-viscosity.sh - Green-Kubo viscosity from viscosity-running.dat
#  usage: lammps-viscosity.sh [-i viscosity-running.dat] [-f 0.5] [-o viscosity]
#  The running Green-Kubo integral should reach a plateau; the value is the mean
#  over the last fraction -f of the run, with its standard deviation as spread.
#  (The original extract-vis-condut averaged after a fixed 100 steps.)
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"
in=viscosity-running.dat; fr=0.5; out=viscosity
while getopts "i:f:o:h" o; do case $o in i) in=$OPTARG ;; f) fr=$OPTARG ;; o) out=$OPTARG ;; *) sed -n '2,8p' "$0"; exit 0 ;; esac; done
[[ -s "$in" ]] || die "no $in"
read -r eta sd n <<< "$(awk '!/^#/ && NF >= 2 { n++; v[n] = $2 } END { s = int(n * (1 - '"$fr"')) + 1; for (i = s; i <= n; i++) { m++; a += v[i]; b += v[i]^2 }
    mu = a/m; printf "%.5g %.3g %d", mu, sqrt(b/m - mu*mu > 0 ? b/m - mu*mu : 0), m }' "$in")"
{ echo "viscosity,std_dev,points"; echo "$eta,$sd,$n"; } > "$out.csv"
info "viscosity = $eta +/- $sd (plateau over $n points); $out.csv"
command -v gnuplot >/dev/null && gnuplot <<EOG 2>/dev/null && info "plot: $out.png"
set terminal pngcairo size 900,650 enhanced font 'Arial,13'
set output '$out.png'
set xlabel 'step'; set ylabel '{/Symbol h} (running Green-Kubo integral)'; set grid lc rgb '#dddddd'; set key off
plot '$in' u 1:2 w l lw 2, $eta lw 1 dt 2
EOG
