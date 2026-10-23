#!/usr/bin/env bash
# =============================================================================
#  lammps-elastic.sh - Young's modulus from the stress-strain curve
#  usage: lammps-elastic.sh [-i stress-strain.dat] [-m 0.02] [-w 51] [-o elastic]
#    -w  running-mean window applied to the stress before fitting (default 51 points)
#  E = slope of a straight line through stress vs strain up to the strain -m
#  (the linear elastic region, default 2 %); the yield point is reported as the
#  maximum stress before it first drops by 5 %.
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"
in=stress-strain.dat; mx=0.02; out=elastic; smooth=51
while getopts "i:m:w:o:h" o; do case $o in i) in=$OPTARG ;; m) mx=$OPTARG ;; w) smooth=$OPTARG ;; o) out=$OPTARG ;; *) sed -n '2,8p' "$0"; exit 0 ;; esac; done
[[ -s "$in" ]] || die "no $in"
res=$(awk -v mx="$mx" -v w="$smooth" -v out="$out" '!/^#/ && NF >= 3 { n++; e[n] = $2; s[n] = $3 }
     END { if (n < 10) { print "ERR fewer than 10 points"; exit }
           # running mean over w points: the raw stress of a small box is very noisy
           for (i = 1; i <= n; i++) { a = i - int(w/2); b = a + w - 1; if (a < 1) { a = 1; b = w } if (b > n) { b = n; a = n - w + 1 }
               t = 0; for (k = a; k <= b; k++) t += s[k]; sm[i] = t / (b - a + 1) }
           for (i = 1; i <= n; i++) if (e[i] <= mx) { m++; sx += e[i]; sy += sm[i]; sxx += e[i]*e[i]; sxy += e[i]*sm[i] }
           if (m < 3) { print "ERR fewer than 3 points below the strain limit"; exit }
           E = (m*sxy - sx*sy)/(m*sxx - sx*sx)
           for (i = 1; i <= n; i++) { if (sm[i] > smax) { smax = sm[i]; emax = e[i] } else if (smax > 0 && sm[i] < 0.95*smax) break }
           for (i = 1; i <= n; i++) print e[i], sm[i] > (out "-smoothed.dat")
           printf "%.5g %.5g %.5g %d", E, smax, emax, m }' "$in")
[[ "$res" == ERR* ]] && die "${res#ERR }"
read -r E smax emax m <<< "$res"
{ echo "Young_modulus,yield_stress,yield_strain,points_in_fit"; echo "$E,$smax,$emax,$m"; } > "$out.csv"
info "E = $E (stress units of the run), yield stress $smax at strain $emax; $out.csv"
command -v gnuplot >/dev/null && gnuplot <<EOG 2>/dev/null && info "plot: $out.png"
set terminal pngcairo size 900,650 enhanced font 'Arial,13'
set output '$out.png'
set xlabel 'strain'; set ylabel 'stress'; set grid lc rgb '#dddddd'; set key top left
plot '$in' u 2:3 w l lc rgb '#cccccc' t 'raw', '$out-smoothed.dat' u 1:2 w l lw 2 lc rgb '#1f5f99' t 'running mean', [0:$mx*2] $E*x t sprintf('E = %.3g', $E) lw 2 dt 2 lc rgb '#b2182b'
EOG
