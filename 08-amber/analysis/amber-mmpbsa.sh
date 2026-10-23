#!/usr/bin/env bash
# =============================================================================
#  amber-mmpbsa.sh - MM-GBSA / MM-PBSA binding free energy with per-residue
#                    decomposition, then a table of the residues that matter
#
#  usage: amber-mmpbsa.sh -d SYSTEM_DIR [-y prod.nc] [-f first:last:interval]
#                         [-g 5|8|2] [-P] [-r "1-50,100-120"] [-t 15]
#    -y  trajectory (default prod.nc)
#    -f  frames (default 1:last:1, "last" = all frames)
#    -g  GB model igb (default 5, OBC II with mbondi2 radii, as in the original work;
#        igb=8 needs mbondi3 radii)
#    -P  also run Poisson-Boltzmann (slower)
#    -r  residues to decompose (default: within 6 A of the ligand in complex-dry.pdb)
#    -t  number of top residues to list (default 15)
#  Uses the dry complex/receptor/ligand topologies written by amber-build.sh, so
#  ante-MMPBSA.py is not needed. Writes mmpbsa/FINAL_RESULTS_MMPBSA.dat,
#  mmpbsa/summary.csv (dG_bind and its terms) and mmpbsa/top-residues.csv.
#  Entropy is not included; relative binding energies across a ligand series
#  are the sound use of these numbers.
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"
d=""; traj=prod.nc; frames="1:last:1"; igb=5; pb=0; resl=""; top=15
while getopts "d:y:f:g:Pr:t:h" o; do
    case $o in d) d=$OPTARG ;; y) traj=$OPTARG ;; f) frames=$OPTARG ;; g) igb=$OPTARG ;; P) pb=1 ;; r) resl=$OPTARG ;; t) top=$OPTARG ;;
               *) sed -n '2,21p' "$0"; exit 0 ;; esac
done
for f in solv.prmtop complex.prmtop receptor.prmtop ligand.prmtop; do [[ -f "$d/$f" ]] || die "$d/$f missing (build with amber-build.sh -p ... -L ...)"; done
[[ -f "$d/$traj" ]] || die "no trajectory $d/$traj"
require MMPBSA.py
IFS=: read -r fs fe fi_ <<< "$frames"
[[ $fe == last ]] && fe=9999999
if [[ -z "$resl" ]]; then
    # residues with any atom within 6 A of the ligand (last residue of the complex is the ligand)
    resl=$(awk '/^(ATOM|HETATM)/ { r = substr($0, 23, 4) + 0; x = substr($0, 31, 8); y = substr($0, 39, 8); z = substr($0, 47, 8); n++; R[n] = r; X[n] = x; Y[n] = y; Z[n] = z; if (r > last) last = r }
        END { for (i = 1; i <= n; i++) if (R[i] == last) { m++; lx[m] = X[i]; ly[m] = Y[i]; lz[m] = Z[i] }
              for (i = 1; i <= n; i++) if (R[i] != last && !(R[i] in near)) for (j = 1; j <= m; j++) if ((X[i]-lx[j])^2 + (Y[i]-ly[j])^2 + (Z[i]-lz[j])^2 < 36) { near[R[i]] = 1; break }
              for (r in near) print r; print last }' "$d/complex-dry.pdb" | sort -n | paste -sd, -)
fi
radii=mbondi2; [[ $igb == 8 ]] && radii=mbondi3
mkdir -p "$d/mmpbsa"
{
    echo "MM-GBSA$( ((pb)) && echo '/PBSA') with per-residue decomposition"
    echo " &general"
    echo "  startframe=$fs, endframe=$fe, interval=${fi_:-1}, keep_files=0, verbose=1,"
    echo "  strip_mask=':WAT:Na+:Cl-:K+',"
    echo " /"
    echo " &gb"; echo "  igb=$igb, saltcon=0.150,"; echo " /"
    ((pb)) && { echo " &pb"; echo "  istrng=0.150, fillratio=4.0, radiopt=0, inp=1,"; echo " /"; }
    echo " &decomp"; echo "  idecomp=2, dec_verbose=0, print_res='$resl',"; echo " /"
} > "$d/mmpbsa/mmpbsa.in"
[[ $igb == 8 ]] && warn "igb=8 needs topologies built with mbondi3 radii (set default PBRadii mbondi3 in leap)"
( cd "$d/mmpbsa" && MMPBSA.py -O -i mmpbsa.in -o FINAL_RESULTS_MMPBSA.dat -do FINAL_DECOMP_MMPBSA.dat \
      -sp ../solv.prmtop -cp ../complex.prmtop -rp ../receptor.prmtop -lp ../ligand.prmtop -y "../$traj" > mmpbsa.log 2>&1 ) ||
    { tail -5 "$d/mmpbsa/mmpbsa.log" >&2; die "MMPBSA.py failed ($d/mmpbsa/mmpbsa.log)"; }
# summary: the "Differences (Complex - Receptor - Ligand)" block of each method
awk 'BEGIN { print "method,term,average_kcal_mol,std_dev,std_err" }
     /GENERALIZED BORN:/ { m = "GB" } /POISSON BOLTZMANN:/ { m = "PB" }
     /Differences \(Complex - Receptor - Ligand\):/ { in_ = 1; next }
     in_ && /^(VDWAALS|EEL|EGB|ESURF|EPB|ENPOLAR|EDISPER|DELTA G gas|DELTA G solv|DELTA TOTAL)/ {
         t = $0; sub(/ +-?[0-9].*$/, "", t); n = split($0, w, " "); printf "%s,%s,%s,%s,%s\n", m, t, w[n-2], w[n-1], w[n] }
     in_ && /^DELTA TOTAL/ { in_ = 0 }' "$d/mmpbsa/FINAL_RESULTS_MMPBSA.dat" > "$d/mmpbsa/summary.csv"
# per-residue total contribution (GB block, DELTAS, "Total Energy Decomposition")
awk -F, '/Generalized Born/ { gb = 1 } /Poisson Boltzmann/ { gb = 0 }
         gb && /DELTAS:/ { de = 1 } gb && de && /Total Energy Decomposition/ { tot = 1; next }
         gb && de && tot && /Sidechain Energy Decomposition/ { exit }
         gb && de && tot && NF > 10 && $1 !~ /Residue/ { print $1 "," $(NF-2) "," $(NF-1) }' "$d/mmpbsa/FINAL_DECOMP_MMPBSA.dat" |
    sort -t, -k2,2g | head -n "$top" | { echo "residue,dG_kcal_mol,std_dev"; cat; } > "$d/mmpbsa/top-residues.csv"
dg=$(awk -F, '$1=="GB" && $2=="DELTA TOTAL"{print $3" +/- "$4}' "$d/mmpbsa/summary.csv")
info "MM-GBSA dG_bind = $dg kcal/mol; tables in $d/mmpbsa/"
