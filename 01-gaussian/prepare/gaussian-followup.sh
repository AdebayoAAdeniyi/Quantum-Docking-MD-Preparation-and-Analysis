#!/usr/bin/env bash
# =============================================================================
#  gaussian-followup.sh - extra jobs that restart from a finished checkpoint
#
#  Reads geometry, basis set (and ECP) from NAME.chk, so nothing has to be
#  re-entered. Run it inside the folder made by gaussian-input.sh (or point
#  -p at it); every job is written next to the checkpoint.
#
#  usage: gaussian-followup.sh -k KIND [options] [-l list | names...]
#    -k temps    thermochemistry at several temperatures from the stored force
#                constants (Freq=ReadFC, no new Hessian)   [-T "283.15 298.15 ..."]
#       solvents TD-DFT in several solvents (solvatochromism) [-S "water,methanol"]
#       nbo      NBO analysis: bond indices and natural resonance theory
#                (needs Gaussian linked to NBO 6/7)
#       wfx      wavefunction file for Multiwfn / AIMAll
#    -m METHOD   must match the method used for NAME.chk (default B3LYP)
#    -D DISP     dispersion (default GD3BJ)
#    -T LIST     temperatures in K (default 283.15 to 333.15 every 5 K)
#    -P ATM      pressure for -k temps (default 1)
#    -S LIST     comma separated solvents (default water,acetonitrile,methanol,dmso,chloroform,cyclohexane)
#    -t N        TD-DFT states (default 20)
#    -x MODEL    solvent model for -k solvents (default SMD)
#    -c FILE     charge/multiplicity table (default chargemult.txt)
#    -p DIR      parent folder holding NAME/NAME.chk (default gaussian-jobs)
#    -n NPROC -M MEM
#  example:
#    gaussian-followup.sh -k temps -l list
#    gaussian-followup.sh -k solvents -S "water,toluene" -l list
#    qdm-submit -s run-temps.sh -l list -p gaussian-jobs/
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"
source "$QDM_ROOT/lib/jobscript.sh"
load_cluster_conf
kind=""; method=B3LYP; disp=GD3BJ; temps=$(seq 283.15 5 333.15 | tr '\n' ' '); press=1
solvents="water,acetonitrile,methanol,dmso,chloroform,cyclohexane"; nstates=20; model=SMD
cmfile=chargemult.txt; parent=gaussian-jobs; nproc=${CORES:-8}; mem=${MEM:-16G}; list=""
while getopts "k:m:D:T:P:S:t:x:c:p:n:M:l:h" o; do
    case $o in
        k) kind=$OPTARG ;; m) method=$OPTARG ;; D) disp=$OPTARG ;; T) temps=$OPTARG ;; P) press=$OPTARG ;;
        S) solvents=$OPTARG ;; t) nstates=$OPTARG ;; x) model=$OPTARG ;; c) cmfile=$OPTARG ;;
        p) parent=$OPTARG ;; n) nproc=$OPTARG ;; M) mem=$OPTARG ;; l) list=$OPTARG ;;
        *) sed -n '2,30p' "$0"; exit 0 ;;
    esac
done
shift $((OPTIND-1))
CORES=$nproc; MEM=$mem; JOB_THREADED=1
[[ "$kind" =~ ^(temps|solvents|nbo|wfx)$ ]] || die "choose -k temps, solvents, nbo or wfx"
names=("$@"); [[ -n "$list" ]] && mapfile -t -O "${#names[@]}" names < <(read_list "$list")
((${#names[@]})) || die "no molecule names given"
dispkw=""; [[ "${disp,,}" != none ]] && dispkw=" EmpiricalDispersion=$disp"
gmem=${mem%B}; gmem=${gmem/%G/GB}; gmem=${gmem/%M/MB}

step() {   # step OLDCHK CHK ROUTE TITLE CHARGE MULT [tail lines]
    printf '%%OldChk=%s\n%%Chk=%s\n%%NProcShared=%s\n%%Mem=%s\n#P %s/ChkBasis%s %s Geom=Check Guess=Read\n\n%s\n\n%s %s\n\n' \
        "$1" "$2" "$nproc" "$gmem" "$method" "$dispkw" "$3" "$4" "$5" "$6"
    [[ -n "${7:-}" ]] && printf '%s\n\n' "$7"
    return 0
}

for name in "${names[@]}"; do
    d="$parent/$name"
    [[ -f "$d/$name.chk" ]] || warn "$name: $d/$name.chk not found yet (input written anyway)"
    q=$(get_charge "$name" "$cmfile"); m=$(get_mult "$name" "$cmfile")
    job="$name-$kind"; mkdir -p "$d"
    {
        i=0
        case $kind in
        temps)
            for T in $temps; do
                ((i++)) && echo "--Link1--"
                step "$name.chk" "$job.chk" "Freq=ReadFC Temperature=$T Pressure=$press" "$name thermochemistry at $T K" "$q" "$m"
            done ;;
        solvents)
            for s in ${solvents//,/ }; do
                ((i++)) && echo "--Link1--"
                step "$name.chk" "$name-td-$s.chk" "TD=(NStates=$nstates) SCRF=($model,Solvent=$s)" "$name TD-DFT in $s" "$q" "$m"
            done ;;
        nbo)
            step "$name.chk" "$job.chk" "Pop=NBORead" "$name NBO bond indices and NRT" "$q" "$m" '$NBO BNDIDX NRT $END' ;;
        wfx)
            step "$name.chk" "$job.chk" "Output=WFX Density=Current" "$name wavefunction" "$q" "$m" "$name.wfx" ;;
        esac
    } > "$d/$job.com"
    write_jobscript "g-$job" "$PWD/$d" GAUSSIAN_SETUP "g16 < $job.com > $job.log" > "$d/run-$kind.sh"
    chmod +x "$d/run-$kind.sh"
    info "$name: $d/$job.com  (run with run-$kind.sh)"
done
