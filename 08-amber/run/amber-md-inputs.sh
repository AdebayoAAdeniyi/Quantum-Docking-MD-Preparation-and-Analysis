#!/usr/bin/env bash
# =============================================================================
#  amber-md-inputs.sh - the standard Amber protocol for a solvated system
#
#   1 min1     solvent and ions relax, solute restrained (10 kcal/mol/A^2)
#   2 min2     everything relaxes
#   3 heat     NVT, 0 -> T K over the first half, solute restrained (5 kcal/mol/A^2)
#   4 density  NPT, restraints 2 kcal/mol/A^2, box adjusts to the right density
#   5 equil    NPT, no restraints
#   6 prod     NPT production, Langevin thermostat, Monte Carlo barostat,
#              SHAKE on H bonds, 2 fs steps, coordinates every -f ps
#
#  usage: amber-md-inputs.sh -d SYSTEM_DIR [-T 300] [-H 100] [-D 100] [-E 1000]
#                            [-P 10] [-f 10] [-m ':LIG|@CA,C,N'] [-c 10] [-g]
#    -d   folder with solv.prmtop / solv.inpcrd (from amber-build.sh)
#    -T   temperature (K)        -H/-D/-E  heat / density / equilibration (ps)
#    -P   production (ns)        -f        frame interval (ps)
#    -m   restraint mask (default: everything that is not water or ions)
#    -c   non-bonded cutoff (A, default 10; PME handles the long range part)
#    -g   use pmemd.cuda (GPU) instead of pmemd.MPI / sander
#  The job script runs the stages in order, stops at the first failure and
#  restarts from the last stage that has a restart file.
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"
source "$QDM_ROOT/lib/jobscript.sh"; load_cluster_conf
d=""; T=300; heat=100; dens=100; equil=1000; prod=10; frame=10; mask='!:WAT,Na+,Cl-,K+'; cut=10; gpu=0
while getopts "d:T:H:D:E:P:f:m:c:gh" o; do
    case $o in d) d=$OPTARG ;; T) T=$OPTARG ;; H) heat=$OPTARG ;; D) dens=$OPTARG ;; E) equil=$OPTARG ;; P) prod=$OPTARG ;;
               f) frame=$OPTARG ;; m) mask=$OPTARG ;; c) cut=$OPTARG ;; g) gpu=1 ;; *) sed -n '2,24p' "$0"; exit 0 ;; esac
done
[[ -f "$d/solv.prmtop" ]] || die "no solv.prmtop in $d"
st() { awk -v ps="$1" 'BEGIN{printf "%d", ps / 0.002}'; }      # ps -> number of 2 fs steps
nw=$(awk -v f="$frame" 'BEGIN{printf "%d", f / 0.002}')
cat > "$d/min1.in" <<EOF
Minimisation 1: solvent and ions, solute restrained
 &cntrl
  imin=1, maxcyc=5000, ncyc=2500, ntb=1, cut=$cut, ntpr=500,
  ntr=1, restraint_wt=10.0, restraintmask='$mask',
 /
EOF
cat > "$d/min2.in" <<EOF
Minimisation 2: whole system
 &cntrl
  imin=1, maxcyc=5000, ncyc=2500, ntb=1, cut=$cut, ntpr=500,
 /
EOF
cat > "$d/heat.in" <<EOF
Heating 0 -> $T K at constant volume, solute restrained
 &cntrl
  imin=0, irest=0, ntx=1, nstlim=$(st "$heat"), dt=0.002,
  ntc=2, ntf=2, cut=$cut, ntb=1,
  ntpr=500, ntwx=$nw, ntwr=5000, ioutfm=1, iwrap=1,
  ntt=3, gamma_ln=2.0, tempi=0.0, temp0=$T, ig=-1,
  ntr=1, restraint_wt=5.0, restraintmask='$mask',
  nmropt=1,
 /
 &wt type='TEMP0', istep1=0, istep2=$(st "$(calc "$heat / 2")"), value1=0.0, value2=$T, /
 &wt type='TEMP0', istep1=$(( $(st "$(calc "$heat / 2")") + 1 )), istep2=$(st "$heat"), value1=$T, value2=$T, /
 &wt type='END' /
EOF
cat > "$d/density.in" <<EOF
Density: NPT with weak restraints
 &cntrl
  imin=0, irest=1, ntx=5, nstlim=$(st "$dens"), dt=0.002,
  ntc=2, ntf=2, cut=$cut, ntb=2, ntp=1, barostat=2, pres0=1.0,
  ntpr=500, ntwx=$nw, ntwr=5000, ioutfm=1, iwrap=1,
  ntt=3, gamma_ln=2.0, temp0=$T, ig=-1,
  ntr=1, restraint_wt=2.0, restraintmask='$mask',
 /
EOF
cat > "$d/equil.in" <<EOF
Equilibration: NPT, no restraints
 &cntrl
  imin=0, irest=1, ntx=5, nstlim=$(st "$equil"), dt=0.002,
  ntc=2, ntf=2, cut=$cut, ntb=2, ntp=1, barostat=2, pres0=1.0,
  ntpr=500, ntwx=$nw, ntwr=5000, ioutfm=1, iwrap=1,
  ntt=3, gamma_ln=2.0, temp0=$T, ig=-1,
 /
EOF
sed -e "s/^Equilibration: NPT, no restraints/Production: NPT, $prod ns/" \
    -e "s/nstlim=$(st "$equil")/nstlim=$(st "$(calc "$prod * 1000")")/" "$d/equil.in" > "$d/prod.in"
cat > "$d/.stages" <<EOF
if (( $gpu )); then EXE=pmemd.cuda
elif command -v pmemd.MPI >/dev/null; then EXE="mpirun -np \$NCORES pmemd.MPI"
elif command -v pmemd >/dev/null; then EXE=pmemd
else EXE=sander; fi
run() {   # run STAGE INPUT_COORDS [restraint reference]
    local s=\$1 c=\$2 ref=\${3:-}
    [[ -s \$s.rst7 ]] && { echo "\$s already done"; return 0; }
    \$EXE -O -i \$s.in -p solv.prmtop -c \$c -o \$s.out -r \$s.rst7 -x \$s.nc -inf \$s.mdinfo \${ref:+-ref \$ref} || return 1
    grep -q "Final Performance Info\\|FINAL RESULTS\\|Total wall time" \$s.out || { echo "\$s did not finish"; return 1; }
}
run min1 solv.inpcrd solv.inpcrd && run min2 min1.rst7 && run heat min2.rst7 min2.rst7 &&
run density heat.rst7 heat.rst7 && run equil density.rst7 && run prod equil.rst7
EOF
write_jobscript "md-$(basename "$d")" "$PWD/$d" AMBER_SETUP "$(cat "$d/.stages")" > "$d/run-md.sh"
rm -f "$d/.stages"; chmod +x "$d/run-md.sh"
info "$d: min1 min2 heat (${heat} ps) density (${dens} ps) equil (${equil} ps) prod (${prod} ns); submit run-md.sh"
