#!/usr/bin/env bash
# =============================================================================
#  plumed-fes.sh - free energy profile from metadynamics hills, with a
#                  convergence check, and the funnel binding free energy
#
#  usage: plumed-fes.sh [-i HILLS] [-k kT] [-n 6] [-u SITE_MAX:UNBOUND_MIN] [-r rcyl] [-o fes]
#    -k   kT in the energy unit of the hills (default 2.494339 kJ/mol = 300 K)
#    -n   number of intermediate profiles for the convergence plot (default 6)
#    -u   funnel analysis: projection limits (nm) of the bound site (z < SITE_MAX)
#         and of the unbound plateau (z > UNBOUND_MIN)
#    -r   funnel cylinder radius (nm, default 0.1)
#  Funnel binding constant (Limongelli, Bonomi and Parrinello, PNAS 2013, 110, 6358):
#    K_b = pi R_cyl^2 * integral_site exp(-[W(z) - W_unbound]/kT) dz
#    dG_b = -kT ln(C0 K_b), C0 = 1/1.66 nm^-3 (1 M)
#  Writes fes.dat, fes-convergence.dat, fes.png and, with -u, binding.csv.
# =============================================================================
source "$(dirname "$0")/../lib/common.sh"
hills=HILLS; kt=2.494339; nint=6; unb=""; rcyl=0.1; out=fes
while getopts "i:k:n:u:r:o:h" o; do case $o in i) hills=$OPTARG ;; k) kt=$OPTARG ;; n) nint=$OPTARG ;; u) unb=$OPTARG ;; r) rcyl=$OPTARG ;; o) out=$OPTARG ;; *) sed -n '2,17p' "$0"; exit 0 ;; esac; done
[[ -s "$hills" ]] || die "no hills file $hills"
require plumed
nh=$(grep -vc '^#' "$hills"); stride=$(( nh / nint )); ((stride > 0)) || stride=1
plumed sum_hills --hills "$hills" --outfile "$out.dat" --mintozero --kt "$kt" > "$out.log" 2>&1 || die "sum_hills failed ($out.log)"
plumed sum_hills --hills "$hills" --outfile "$out-part-" --mintozero --kt "$kt" --stride "$stride" >> "$out.log" 2>&1 || true
: > "$out-convergence.dat"
for f in $(ls "$out-part-"*.dat 2>/dev/null | sort -V); do grep -v '^#' "$f" | awk '{print $1, $2}' >> "$out-convergence.dat"; printf '\n\n' >> "$out-convergence.dat"; done
rm -f "$out-part-"*.dat
if [[ -n "$unb" ]]; then
    IFS=: read -r zs zu <<< "$unb"
    awk -v zs="$zs" -v zu="$zu" -v kt="$kt" -v r="$rcyl" '!/^#/ && NF >= 2 { n++; z[n] = $1; w[n] = $2 }
        END { for (i = 1; i <= n; i++) if (z[i] > zu) { s += w[i]; m++ }
              if (!m) { print "no unbound points beyond " zu > "/dev/stderr"; exit 1 }
              wu = s / m; dz = z[2] - z[1]
              for (i = 1; i <= n; i++) if (z[i] < zs) I += exp(-(w[i] - wu) / kt) * dz
              pi = 3.141592653589793; Kb = pi * r * r * I; C0 = 1 / 1.6605
              printf "W_unbound_kJ_mol,K_b_nm3,dG_binding_kJ_mol,dG_binding_kcal_mol\n%.3f,%.5g,%.2f,%.2f\n", wu, Kb, -kt * log(C0 * Kb), -kt * log(C0 * Kb) / 4.184 }' "$out.dat" > "binding.csv" &&
        info "funnel binding free energy: $(tail -1 binding.csv | cut -d, -f3) kJ/mol (binding.csv)"
fi
if command -v gnuplot >/dev/null; then gnuplot <<EOG 2>/dev/null && info "plot: $out.png"
set terminal pngcairo size 900,650 enhanced font 'Arial,13'
set output '$out.png'
set xlabel 'collective variable'; set ylabel 'free energy'; set grid lc rgb '#dddddd'; set key top right
plot for [i=0:$((nint - 1))] '$out-convergence.dat' index i u 1:2 w l lc rgb '#bbbbbb' notitle, '$out.dat' u 1:2 w l lw 3 lc rgb '#1f5f99' t 'final'
EOG
fi
info "FES: $out.dat ($nh hills); convergence: $out-convergence.dat"
