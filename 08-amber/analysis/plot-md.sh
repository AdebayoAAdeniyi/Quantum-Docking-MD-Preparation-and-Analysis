#!/usr/bin/env bash
# =============================================================================
#  plot-md.sh - PNG plots of the usual MD time series in an analysis folder
#  usage: plot-md.sh ANALYSIS_DIR
#  Recognises (cpptraj or gmx output, two columns, # comments allowed):
#   energy.dat rmsd-backbone.dat rmsd-ligand.dat rmsf-ca.dat radgyr.dat
#   hbonds-solute.dat (number of H-bonds per frame), sasa.dat
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"; require gnuplot
a=${1:?analysis folder}
plot1() {   # plot1 FILE XLABEL YLABEL TITLE [COL]
    [[ -s "$a/$1" ]] || return 0
    gnuplot <<EOG 2>/dev/null
set terminal pngcairo size 900,600 enhanced font 'Arial,12'
set output '$a/${1%.*}.png'
set xlabel '$2'; set ylabel '$3'; set title '$4'; set grid lc rgb '#dddddd'; set key off
set datafile commentschars '#@&'
plot '$a/$1' u 1:${5:-2} w l lw 1.5 lc rgb '#1f5f99'
EOG
    info "plot: $a/${1%.*}.png"
}
plot1 rmsd-backbone.dat 'frame' 'RMSD (Å)' 'Protein backbone RMSD'
plot1 rmsd-ligand.dat   'frame' 'RMSD (Å)' 'Ligand RMSD (heavy atoms, protein frame)'
plot1 rmsf-ca.dat       'residue' 'RMSF (Å)' 'C{/Symbol a} fluctuations'
plot1 radgyr.dat        'frame' 'R_g (Å)' 'Radius of gyration'
plot1 hbonds-solute.dat 'frame' 'number of H-bonds' 'Solute hydrogen bonds'
plot1 sasa.dat          'time' 'SASA' 'Solvent accessible surface area'
if [[ -s "$a/energy.dat" ]]; then
    gnuplot <<EOG 2>/dev/null
set terminal pngcairo size 1000,800 enhanced font 'Arial,12'
set output '$a/energy.png'
set multiplot layout 3,1
set grid lc rgb '#dddddd'; set key off; set xlabel 'time (ps)'
set ylabel 'E_{tot} (kcal/mol)'; plot '$a/energy.dat' u 2:3 w l lc rgb '#1f5f99'
set ylabel 'T (K)'; plot '$a/energy.dat' u 2:5 w l lc rgb '#b2182b'
set ylabel 'density (g/cm^3)'; plot '$a/energy.dat' u 2:(\$6 > 0 ? \$6 : 1/0) w l lc rgb '#1b7837'
unset multiplot
EOG
    info "plot: $a/energy.png"
fi
