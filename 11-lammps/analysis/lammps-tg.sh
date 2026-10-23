#!/usr/bin/env bash
# =============================================================================
#  lammps-tg.sh - glass transition temperature from a cooling run
#  usage: lammps-tg.sh [-i tg-density.dat] [-o tg]
#  Input: "T <T> <density>" lines written by glass-transition.in.
#  Two straight lines are fitted to specific volume (1/density) against T, one
#  to the low-T (glassy) and one to the high-T (rubbery/melt) points; every split
#  with at least three points per branch is tried and the split with the lowest
#  total squared error is kept. Tg is where the two lines cross. Also reports the
#  thermal expansion coefficient of each branch, alpha = (dV/dT)/V.
#  The original plot-glass did the same fit with a fixed, hand-chosen split.
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"
in=tg-density.dat; out=tg
while getopts "i:o:h" o; do case $o in i) in=$OPTARG ;; o) out=$OPTARG ;; *) sed -n '2,13p' "$0"; exit 0 ;; esac; done
[[ -s "$in" ]] || die "no $in"
sort -g -k1,1 "$in" | awk '$3 > 0 {print $1, 1/$3}' > "$out-volume.dat"
res=$(awk '{ n++; T[n] = $1; V[n] = $2 }
    function fit(a, b,   i, m, sx, sy, sxx, sxy, e) { m = sx = sy = sxx = sxy = 0
        for (i = a; i <= b; i++) { m++; sx += T[i]; sy += V[i]; sxx += T[i]^2; sxy += T[i]*V[i] }
        S = (m*sxy - sx*sy) / (m*sxx - sx*sx); I = (sy - S*sx)/m; e = 0
        for (i = a; i <= b; i++) e += (V[i] - S*T[i] - I)^2; return e }
    END { if (n < 6) { print "ERR need at least 6 temperatures"; exit }
          best = -1; for (k = 3; k <= n - 3; k++) { e = fit(1, k) + fit(k + 1, n); if (best < 0 || e < best) { best = e; kb = k } }
          fit(1, kb); s1 = S; i1 = I; fit(kb + 1, n); s2 = S; i2 = I
          tg = (i1 - i2) / (s2 - s1); vg = s1 * tg + i1
          printf "%.4f %.6g %.6g %.6g %.6g %.4f %.4f\n", tg, s1, i1, s2, i2, s1 / vg, s2 / vg }' "$out-volume.dat")
[[ "$res" == ERR* ]] && die "${res#ERR }"
read -r tg s1 i1 s2 i2 a1 a2 <<< "$res"
{ echo "Tg,slope_glass,intercept_glass,slope_melt,intercept_melt,alpha_glass,alpha_melt"; echo "$tg,$s1,$i1,$s2,$i2,$a1,$a2" | tr ' ' ','; } > "$out.csv"
info "Tg = $tg (alpha glass $a1, melt $a2 per T unit); $out.csv"
if command -v gnuplot >/dev/null; then gnuplot <<EOG 2>/dev/null
set terminal pngcairo size 900,650 enhanced font 'Arial,13'
set output '$out.png'
set xlabel 'T'; set ylabel 'specific volume (1/{/Symbol r})'; set grid lc rgb '#dddddd'; set key top left
set arrow from $tg, graph 0 to $tg, graph 1 nohead dt 2 lc rgb 'grey'
set label sprintf('T_g = %.3f', $tg) at $tg, graph 0.9 offset 1,0
plot '$out-volume.dat' u 1:2 w p pt 7 ps 1.3 t 'simulation', $s1*x+$i1 t 'glass fit' lw 2, $s2*x+$i2 t 'melt fit' lw 2
EOG
info "plot: $out.png"; fi
