#!/usr/bin/env bash
# =============================================================================
#  xtb-run.sh - fast GFN2-xTB screening of a whole set of molecules
#
#  Semi-empirical tight binding (Bannwarth, Ehlert and Grimme, J. Chem. Theory
#  Comput. 2019, 15, 1652) takes seconds per molecule, so it is used here to
#  pre-optimise structures before DFT, to prepare docking ligands, and to rank
#  a series by gap, IP, EA and reactivity indices before spending DFT time.
#  It replaces the old "quick PM3 pre-optimisation through Gaussian under wine".
#
#  usage: xtb-run.sh [-j STEPS] [-s SOLVENT] [-c chargemult.txt] [-o OUTDIR]
#                    [-q] [-n THREADS] [-l list -d dir | files...]
#    -j  comma separated (default ohess,vipea,vfukui):
#          opt     optimisation only            ohess  optimisation + Hessian,
#          vipea   vertical IP and EA                  thermochemistry, IR
#          vfukui  condensed Fukui indices
#    -s  ALPB implicit solvent (water, methanol, acetonitrile, dmso, thf ...)
#    -q  only write the job scripts (run later with qdm-submit); by default
#        every molecule is run immediately on this machine
#  Later steps use the optimised geometry (NAME-xtbopt.xyz).
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"
source "$QDM_ROOT/lib/jobscript.sh"
load_cluster_conf
steps=ohess,vipea,vfukui; solv=""; cmfile=chargemult.txt; outdir=xtb-jobs; queue=0; nthr=${CORES:-4}; list=""; indir=.
while getopts "j:s:c:o:qn:l:d:h" o; do
    case $o in j) steps=$OPTARG ;; s) solv=$OPTARG ;; c) cmfile=$OPTARG ;; o) outdir=$OPTARG ;; q) queue=1 ;;
               n) nthr=$OPTARG ;; l) list=$OPTARG ;; d) indir=$OPTARG ;; *) sed -n '2,22p' "$0"; exit 0 ;; esac
done
shift $((OPTIND-1))
((queue)) || require xtb
CORES=$nthr; JOB_THREADED=1
mapfile -t files < <(gather_structures "$list" "$indir" "$@")
((${#files[@]})) || die "no structures"
sol=""; [[ -n "$solv" ]] && sol="--alpb $solv"
mkdir -p "$outdir"
for f in "${files[@]}"; do
    n=$(strip_ext "$f"); d="$outdir/$n"; mkdir -p "$d"
    q=$(get_charge "$n" "$cmfile"); m=$(get_mult "$n" "$cmfile"); uhf=$((m - 1))
    to_format "$f" "$d/$n.xyz"; check_multiplicity "$d/$n.xyz" "$q" "$m"
    cmds=("export OMP_NUM_THREADS=$nthr OMP_STACKSIZE=1G" "ulimit -s unlimited 2>/dev/null || true")
    g="$n.xyz"
    for s in ${steps//,/ }; do
        case $s in
            opt|ohess) cmds+=("xtb $g --$s --chrg $q --uhf $uhf $sol > $n-$s.out 2>&1"
                              "cp xtbopt.xyz $n-xtbopt.xyz"
                              "[ -f vibspectrum ] && mv vibspectrum $n-vibspectrum || true")
                       g="$n-xtbopt.xyz" ;;
            vipea|vfukui) cmds+=("xtb $g --$s --chrg $q --uhf $uhf $sol > $n-$s.out 2>&1") ;;
            *) die "unknown step $s" ;;
        esac
    done
    write_jobscript "xtb-$n" "$PWD/$d" XTB_SETUP "${cmds[@]}" > "$d/run.sh"; chmod +x "$d/run.sh"
    if ((queue)); then info "$n: $d/run.sh written"
    else
        if (cd "$d" && SCHEDULER=local bash run.sh > run.log 2>&1); then info "$n: done ($steps)"
        else warn "$n: xtb failed, see $d/run.log"; fi
    fi
done
