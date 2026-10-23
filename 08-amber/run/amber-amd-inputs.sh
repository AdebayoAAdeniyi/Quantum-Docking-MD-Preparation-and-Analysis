#!/usr/bin/env bash
# =============================================================================
#  amber-amd-inputs.sh - accelerated MD (dual boost) parameters from a
#                        conventional MD run, and the aMD production input
#
#  Rule of thumb from the Amber aMD tutorial (Pierce, Salomon-Ferrer, de Oliveira,
#  McCammon and Walker, J. Chem. Theory Comput. 2012, 8, 2997):
#    EthreshP = <EPtot> + 0.16 * N_atoms        alphaP = 0.16 * N_atoms
#    EthreshD = <DIHED> + 4 * N_solute_res      alphaD = 0.2 * 4 * N_solute_res
#  averages taken from the cMD output (default prod.out).
#
#  usage: amber-amd-inputs.sh -d SYSTEM_DIR [-o prod.out] [-n 500] [-T 300]
#    -n   aMD length in ns (default 500 ns split in 100 ns chunks)
#  Writes amd.in and run-amd.sh (pmemd/pmemd.cuda, -amdlog for reweighting).
#  Reweighting of aMD free energy surfaces: PyReweighting (Miao et al., J. Chem.
#  Theory Comput. 2014, 10, 2677), which is not redistributed here.
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"
source "$QDM_ROOT/lib/jobscript.sh"; load_cluster_conf
d=""; mdout=prod.out; ns=500; T=300
while getopts "d:o:n:T:h" o; do case $o in d) d=$OPTARG ;; o) mdout=$OPTARG ;; n) ns=$OPTARG ;; T) T=$OPTARG ;; *) sed -n '2,17p' "$0"; exit 0 ;; esac; done
[[ -f "$d/$mdout" && -f "$d/solv.prmtop" ]] || die "need $d/$mdout and $d/solv.prmtop"
natom=$(awk '/%FLAG POINTERS/{getline; getline; print $1; exit}' "$d/solv.prmtop")
nres=$(awk '/%FLAG RESIDUE_LABEL/{f=1; getline; next} /%FLAG/{f=0} f' "$d/solv.prmtop" | tr -s ' ' '\n' | grep -cvxE 'WAT|HOH|Na\+|Cl-|K\+|')
read -r ep dih n <<< "$(awk '/A V E R A G E S/{exit} / EPtot      =/{ep += $NF; n++} / DIHED      =/{dh += $NF}
                             END { if (n) printf "%.2f %.2f %d", ep/n, dh/n, n }' "$d/$mdout")"
[[ -n "$n" ]] || die "no energies in $d/$mdout"
read -r ethp alpp ethd alpd <<< "$(awk -v ep="$ep" -v dh="$dih" -v na="$natom" -v nr="$nres" \
    'BEGIN { printf "%.1f %.1f %.1f %.1f", ep + 0.16 * na, 0.16 * na, dh + 4 * nr, 0.2 * 4 * nr }')"
cat > "$d/amd.in" <<EOF
Accelerated MD, dual boost (iamd=3), 100 ns per chunk
 &cntrl
  imin=0, irest=1, ntx=5, nstlim=50000000, dt=0.002,
  ntc=2, ntf=2, cut=10.0, ntb=2, ntp=1, barostat=2,
  ntpr=5000, ntwx=5000, ntwr=50000, ioutfm=1, iwrap=1,
  ntt=3, gamma_ln=2.0, temp0=$T, ig=-1,
  iamd=3, ethreshd=$ethd, alphad=$alpd, ethreshp=$ethp, alphap=$alpp,
 /
EOF
chunks=$(( (ns + 99) / 100 ))
cmds=('EXE=pmemd; command -v pmemd.cuda >/dev/null && EXE=pmemd.cuda' 'prev=prod.rst7')
for i in $(seq 1 $chunks); do cmds+=("[[ -s amd-$i.rst7 ]] || \$EXE -O -i amd.in -p solv.prmtop -c \$prev -o amd-$i.out -r amd-$i.rst7 -x amd-$i.nc -amdlog amd-$i.log" "prev=amd-$i.rst7"); done
write_jobscript "amd-$(basename "$d")" "$PWD/$d" AMBER_SETUP "${cmds[@]}" > "$d/run-amd.sh"; chmod +x "$d/run-amd.sh"
info "$d: <EPtot> = $ep, <DIHED> = $dih over $n frames; N_atoms = $natom, solute residues = $nres"
info "EthreshP = $ethp, alphaP = $alpp, EthreshD = $ethd, alphaD = $alpd -> $d/amd.in, run-amd.sh ($chunks x 100 ns)"
