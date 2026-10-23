#!/usr/bin/env bash
# =============================================================================
#  gmx-analysis.sh - standard GROMACS trajectory analysis in one call
#
#  usage: gmx-analysis.sh -d SYSTEM_DIR [-s md] [-b 0] [-k 0.1] [-I]
#    -s  run name (default md -> md.tpr, md.xtc, md.edr)
#    -b  first time to analyse (ps)
#    -k  RMSD cut-off for gromos clustering of the backbone (nm, default 0.1)
#    -I  protein-ligand interaction energy (Coulomb-SR + LJ-SR) by rerun
#  Output in SYSTEM_DIR/analysis (xvg converted to two-column .dat files):
#    energy (potential, temperature, pressure, density), rmsd-backbone, rmsd-ligand
#    (fitted on the backbone), rmsf-ca (per residue), radgyr, hbonds-solute
#    (protein-ligand H-bonds when a ligand exists), sasa, clusters, and the free
#    energy landscape fes-rmsd-rg.dat = -kT ln P(RMSD, Rg) with its contour plot.
#  The old run-gromac-analysis-* scripts (one copy per cluster) looked up energy
#  term numbers by hand; terms are selected by name here.
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"
d=""; s=md; b=0; cut=0.1; inter=0
while getopts "d:s:b:k:Ih" o; do case $o in d) d=$OPTARG ;; s) s=$OPTARG ;; b) b=$OPTARG ;; k) cut=$OPTARG ;; I) inter=1 ;; *) sed -n '2,18p' "$0"; exit 0 ;; esac; done
[[ -f "$d/$s.tpr" && -f "$d/$s.xtc" ]] || die "need $d/$s.tpr and $d/$s.xtc"
GMX=${GMX:-gmx}; a="$d/analysis"; mkdir -p "$a"
ndx="$d/index.ndx"; [[ -f "$ndx" ]] || die "no index.ndx in $d"
has() { grep -q "^\[ $1 \]" "$ndx"; }
lig=""; has LIG && lig=LIG
cplx=Protein; [[ -n "$lig" ]] && cplx=Protein_LIG
xvg2dat() { grep -v '^[#@]' "$1" > "$2"; }
q() { "$GMX" "$@" >> "$a/gmx.log" 2>&1; }
: > "$a/gmx.log"
# 1. whole molecules, solute centred, no jumps
printf '%s\n%s\n' "$cplx" System | q trjconv -s "$d/$s.tpr" -f "$d/$s.xtc" -n "$ndx" -o "$a/center.xtc" -pbc mol -center -ur compact -b "$b" ||
    die "trjconv failed (see $a/gmx.log)"
T="$a/center.xtc"
# 2. energies by name
for term in Potential Temperature Pressure Density; do
    echo "$term" | q energy -f "$d/$s.edr" -o "$a/e-$term.xvg" -b "$b" && xvg2dat "$a/e-$term.xvg" "$a/energy-$term.dat"
done
paste "$a/energy-Potential.dat" "$a/energy-Temperature.dat" "$a/energy-Density.dat" 2>/dev/null |
    awk '{print "md", $1, $2, $2, $4, $6}' > "$a/energy.dat"
# 3. structure
printf 'Backbone\nBackbone\n' | q rms -s "$d/$s.tpr" -f "$T" -n "$ndx" -o "$a/rmsd.xvg" -tu ns && xvg2dat "$a/rmsd.xvg" "$a/rmsd-backbone.dat"
[[ -n "$lig" ]] && printf 'Backbone\nLIG\n' | q rms -s "$d/$s.tpr" -f "$T" -n "$ndx" -o "$a/rmsd-lig.xvg" -tu ns && xvg2dat "$a/rmsd-lig.xvg" "$a/rmsd-ligand.dat"
echo C-alpha | q rmsf -s "$d/$s.tpr" -f "$T" -n "$ndx" -o "$a/rmsf.xvg" -res && xvg2dat "$a/rmsf.xvg" "$a/rmsf-ca.dat"
echo Protein | q gyrate -s "$d/$s.tpr" -f "$T" -n "$ndx" -o "$a/gyrate.xvg" && xvg2dat "$a/gyrate.xvg" "$a/radgyr.dat"
echo Protein | q sasa -s "$d/$s.tpr" -f "$T" -n "$ndx" -o "$a/sasa.xvg" && xvg2dat "$a/sasa.xvg" "$a/sasa.dat"
# 4. hydrogen bonds (GROMACS 2023+ syntax, falling back to the legacy tool)
if [[ -n "$lig" ]]; then
    q hbond -s "$d/$s.tpr" -f "$T" -n "$ndx" -r Protein -t LIG -num "$a/hbnum.xvg" ||
        printf 'Protein\nLIG\n' | q hbond-legacy -s "$d/$s.tpr" -f "$T" -n "$ndx" -num "$a/hbnum.xvg"
else
    q hbond -s "$d/$s.tpr" -f "$T" -n "$ndx" -r Protein -t Protein -num "$a/hbnum.xvg" ||
        printf 'Protein\nProtein\n' | q hbond-legacy -s "$d/$s.tpr" -f "$T" -n "$ndx" -num "$a/hbnum.xvg"
fi
[[ -f "$a/hbnum.xvg" ]] && xvg2dat "$a/hbnum.xvg" "$a/hbonds-solute.dat"
# 5. clustering of the backbone
printf 'Backbone\n%s\n' "$cplx" | q cluster -s "$d/$s.tpr" -f "$T" -n "$ndx" -method gromos -cutoff "$cut" -g "$a/cluster.log" -cl "$a/clusters.pdb" -dist "$a/rmsd-dist.xvg" || true
ncl=$(grep -m1 "Found .* clusters" "$a/cluster.log" 2>/dev/null | awk '{print $2}')
# 6. free energy landscape from RMSD and Rg: G = -kT ln(P/Pmax), kT at 300 K in kJ/mol
if [[ -s "$a/rmsd-backbone.dat" && -s "$a/radgyr.dat" ]]; then
    paste <(awk '{print $2}' "$a/rmsd-backbone.dat") <(awk '{print $2}' "$a/radgyr.dat") | awk -v nb=30 '
        NF == 2 { n++; x[n] = $1; y[n] = $2; if (n == 1 || $1 < xl) xl = $1; if (n == 1 || $1 > xh) xh = $1; if (n == 1 || $2 < yl) yl = $2; if (n == 1 || $2 > yh) yh = $2 }
        END { dx = (xh - xl) / nb + 1e-9; dy = (yh - yl) / nb + 1e-9
              for (i = 1; i <= n; i++) { h[int((x[i] - xl) / dx), int((y[i] - yl) / dy)]++ }
              for (k in h) if (h[k] > m) m = h[k]
              for (i = 0; i < nb; i++) { for (j = 0; j < nb; j++) { c = h[i, j]
                  printf "%.4f %.4f %s\n", xl + (i + 0.5) * dx, yl + (j + 0.5) * dy, (c > 0 ? sprintf("%.3f", -2.494339 * log(c / m)) : "NaN") } print "" } }' > "$a/fes-rmsd-rg.dat"
fi
# 7. protein-ligand interaction energy by rerun with energy groups (CPU)
if ((inter)) && [[ -n "$lig" ]]; then
    sed -e '/^energygrps/d' "$d/$s.mdp" > "$a/rerun.mdp"; echo "energygrps = Protein LIG" >> "$a/rerun.mdp"
    q grompp -f "$a/rerun.mdp" -c "$d/$s.gro" -p "$d/topol.top" -n "$ndx" -o "$a/rerun.tpr" -maxwarn 2 &&
    q mdrun -s "$a/rerun.tpr" -rerun "$d/$s.xtc" -deffnm "$a/rerun" -nb cpu &&
    printf 'Coul-SR:Protein-LIG\nLJ-SR:Protein-LIG\n' | q energy -f "$a/rerun.edr" -o "$a/interaction.xvg" &&
    xvg2dat "$a/interaction.xvg" "$a/interaction-energy.dat" &&
    awk '{c += $2; l += $3; n++} END {printf "mean Coulomb-SR %.1f, LJ-SR %.1f, total %.1f kJ/mol over %d frames\n", c/n, l/n, (c+l)/n, n}' "$a/interaction-energy.dat" > "$a/interaction-summary.txt"
fi
# plots
"$QDM_ROOT/08-amber/analysis/plot-md.sh" "$a" || true
if [[ -s "$a/fes-rmsd-rg.dat" ]] && command -v gnuplot >/dev/null; then
    gnuplot <<EOG 2>/dev/null
set terminal pngcairo size 900,700 enhanced font 'Arial,12'
set output '$a/fes-rmsd-rg.png'
set xlabel 'backbone RMSD (nm)'; set ylabel 'R_g (nm)'; set cblabel 'G (kJ/mol)'
set view map; set palette rgb 33,13,10; set datafile missing 'NaN'
splot '$a/fes-rmsd-rg.dat' u 1:2:3 w pm3d notitle
EOG
fi
info "analysis in $a (${ncl:-?} clusters at $cut nm; ligand: ${lig:-none})"
