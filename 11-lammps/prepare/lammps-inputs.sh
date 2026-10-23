#!/usr/bin/env bash
# =============================================================================
#  lammps-inputs.sh - polymer and materials property runs from one system file
#
#  usage: lammps-inputs.sh -s system.in -D system.data [-u real|lj] [-T T] [-P P]
#                          [-j equilibrate,tg,elastic,viscosity,rdf] [-H Thigh -L Tlow -n nT]
#                          [-e erate] [-x STEPS_SCALE] [-o OUTDIR]
#    system.in   units, styles, "read_data ${data}" and all force field coefficients
#                (the original templates repeated these in every input; now they live
#                in one file and every template includes it)
#    -u          unit system for default timestep, damping and conversion factors
#                (real: fs, atm -> MPa and Pa s; lj: reduced units)
#    -x          multiply all run lengths (default 1; use 0.1 for a quick test)
#  The stages chain through equilibrated.data: equilibrate first, then the others.
#  Writes OUTDIR/params.in, the five inputs and run.sh.
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"; source "$QDM_ROOT/lib/jobscript.sh"; load_cluster_conf
sys=""; data=""; units=real; T=""; P=""; jobs=equilibrate,tg,elastic,viscosity,rdf; th=""; tl=""; nT=11; erate=""; x=1; out=lammps-run
while getopts "s:D:u:T:P:j:H:L:n:e:x:o:h" o; do
    case $o in s) sys=$OPTARG ;; D) data=$OPTARG ;; u) units=$OPTARG ;; T) T=$OPTARG ;; P) P=$OPTARG ;; j) jobs=$OPTARG ;;
               H) th=$OPTARG ;; L) tl=$OPTARG ;; n) nT=$OPTARG ;; e) erate=$OPTARG ;; x) x=$OPTARG ;; o) out=$OPTARG ;; *) sed -n '2,18p' "$0"; exit 0 ;; esac
done
[[ -f "$sys" && -f "$data" ]] || die "need -s system.in and -D data file"
mkdir -p "$out"; cp "$sys" "$out/system.in"; cp "$data" "$out/initial.data"
cp "$QDM_ROOT"/11-lammps/templates/*.in "$out/"
n() { awk -v a="$1" -v x="$x" 'BEGIN{v=int(a*x); if (v<100) v=100; print v}'; }
case $units in
    real) : "${T:=300}" "${P:=1.0}" "${th:=600}" "${tl:=200}" "${erate:=1e-6}"
          dt=1.0; td=100.0; pd=1000.0; pconv=0.101325; kB=1.3806504e-23; gk="101325.0*101325.0*1.0e-15*1.0e-30"
          punit="MPa"; vunit="Pa s"; base=1000000 ;;
    lj)   : "${T:=1.0}" "${P:=0.0}" "${th:=1.0}" "${tl:=0.2}" "${erate:=1e-3}"
          dt=0.005; td=0.5; pd=5.0; pconv=1.0; kB=1.0; gk=1.0; punit="eps/sigma^3"; vunit="reduced"; base=200000 ;;
    *) die "units must be real or lj" ;;
esac
cat > "$out/params.in" <<EOF2
# written by lammps-inputs.sh ($units units); any value can be overridden with -var NAME VALUE
variable data   index initial.data
variable T      index $T
variable P      index $P
variable dt     index $dt
variable tdamp  index $td
variable pdamp  index $pd
variable thermo index 1000
variable nheat  index $(n $((base / 4)))
variable nsteps index $(n $base)
variable Thigh  index $th
variable Tlow   index $tl
variable nT     index $nT
variable nequil index $(n $((base / 10)))
variable nprod  index $(n $((base / 10)))
variable erate  index $erate
variable pconv  index $pconv
variable kB     index $kB
variable gkconv index $gk
variable s      index 5
variable p      index 200
variable d      index 1000
EOF2
LMP=${LMP:-lmp}
cmds=("LMP=\"$LMP\"")
for j in ${jobs//,/ }; do
    case $j in
        equilibrate) cmds+=("[[ -s equilibrated.data ]] || \$LMP -in equilibrate.in -log equilibrate.log") ;;
        tg)          cmds+=("rm -f tg-density.dat; \$LMP -in glass-transition.in -var data equilibrated.data -log tg.log") ;;
        elastic)     cmds+=("\$LMP -in elastic.in -var data equilibrated.data -log elastic.log") ;;
        viscosity)   cmds+=("\$LMP -in viscosity.in -var data equilibrated.data -log viscosity.log") ;;
        rdf)         cmds+=("\$LMP -in rdf-msd.in -var data equilibrated.data -log rdf-msd.log") ;;
        *) die "unknown job $j" ;;
    esac
done
write_jobscript "lmp-$(basename "$out")" "$(cd "$out" && pwd)" LAMMPS_SETUP "${cmds[@]}" > "$out/run.sh"; chmod +x "$out/run.sh"
info "$out: $jobs ($units units, T = $T, P = $P; stress in $punit, viscosity in $vunit)"
