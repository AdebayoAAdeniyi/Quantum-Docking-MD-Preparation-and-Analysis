#!/usr/bin/env bash
# =============================================================================
#  gaussian-input.sh - build multi-step Gaussian 09/16 jobs for a set of molecules
#
#  One input per molecule, made of --Link1-- steps that share one checkpoint,
#  so a single submission gives the optimised structure and all its properties.
#  Any structure format Open Babel reads is accepted (xyz, pdb, mol2, sdf, smi,
#  a previous Gaussian .log ...). Atoms are reordered with metals first and
#  hydrogens last, and an ECP basis is added automatically for heavy atoms.
#
#  usage: gaussian-input.sh [options] [-l list -d dir | files...]
#    -j STEPS   comma separated steps (default opt,pop,chelpg,td,nmr):
#                 opt     optimisation + frequencies (IR, thermochemistry)
#                 optraman  same, with Raman activities
#                 sp      single point energy
#                 pop     NPA + Hirshfeld/CM5 charges, full orbital print
#                 chelpg  electrostatic potential fitted charges (MEP)
#                 td      TD-DFT excited states (UV-Vis)        [-t nstates]
#                 nmr     GIAO shielding (+ spin-spin with -J)
#                 polar   static polarisability alpha and hyperpolarisability beta
#                 epr     hyperfine couplings, g tensor (open shell only)
#                 wfx     wavefunction file for Multiwfn / AIMAll
#               Separate inputs for conceptual DFT are made with -F (see below).
#    -m METHOD  functional or method          (default B3LYP)
#    -D DISP    empirical dispersion: GD3BJ, GD3, none (default GD3BJ)
#    -b BASIS   basis for light atoms          (default 6-311+G(d,p))
#    -e ECP     basis + ECP for heavy atoms    (default LanL2DZ; SDD, def2TZVP ...)
#    -s SOLV    solvent for property steps as MODEL:SOLVENT, e.g. smd:water,
#               cpcm:acetonitrile              (default: gas phase)
#    -O         also optimise in the solvent (default: optimise in gas phase)
#    -t N       number of TD-DFT states        (default 20)
#    -J         add spin-spin couplings to the NMR step
#    -F         also write NAME-fukui.com: N, N+1, N-1 single points at the
#               optimised geometry for Fukui functions, IP, EA, hardness ...
#    -c FILE    charge/multiplicity table "name charge mult" (default chargemult.txt)
#    -n NPROC   processors (default CORES from config/cluster.conf)
#    -M MEM     memory     (default MEM from config/cluster.conf)
#    -o DIR     output directory (default gaussian-jobs)
#  example:
#    gaussian-input.sh -l list -d structures -j opt,pop,td -s smd:acetonitrile -F
#  then:
#    qdm-submit -l list -p gaussian-jobs/
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"
source "$QDM_ROOT/lib/jobscript.sh"
load_cluster_conf

steps=opt,pop,chelpg,td,nmr; method=B3LYP; disp=GD3BJ; basis="6-311+G(d,p)"; ecp=LanL2DZ
solv=""; optsolv=0; nstates=20; spinspin=0; fukui=0; cmfile=chargemult.txt
nproc=${CORES:-8}; mem=${MEM:-16G}; outdir=gaussian-jobs; list=""; indir=.
while getopts "j:m:D:b:e:s:Ot:JFc:n:M:o:l:d:h" o; do
    case $o in
        j) steps=$OPTARG ;; m) method=$OPTARG ;; D) disp=$OPTARG ;; b) basis=$OPTARG ;;
        e) ecp=$OPTARG ;; s) solv=$OPTARG ;; O) optsolv=1 ;; t) nstates=$OPTARG ;;
        J) spinspin=1 ;; F) fukui=1 ;; c) cmfile=$OPTARG ;; n) nproc=$OPTARG ;;
        M) mem=$OPTARG ;; o) outdir=$OPTARG ;; l) list=$OPTARG ;; d) indir=$OPTARG ;;
        *) sed -n '2,43p' "$0"; exit 0 ;;
    esac
done
CORES=$nproc; MEM=$mem; JOB_THREADED=1         # the job script header follows -n / -M
shift $((OPTIND-1))
mapfile -t files < <(gather_structures "$list" "$indir" "$@")
((${#files[@]})) || die "no structures given (use -l list -d dir or file names)"

# ---------- route pieces ----------
dispkw=""; [[ "${disp,,}" != none ]] && dispkw=" EmpiricalDispersion=$disp"
scrf=""
if [[ -n "$solv" ]]; then
    model=${solv%%:*}; solvent=${solv#*:}
    case "${model,,}" in smd|cpcm|pcm|iefpcm) ;; *) die "solvent model must be smd, cpcm or pcm (got $model)" ;; esac
    scrf=" SCRF=(${model^^},Solvent=$solvent)"
fi
optscrf=""; ((optsolv)) && optscrf=$scrf
# G16 accepts "Mem=16GB"; Slurm style "16G" is converted
gmem=${mem%B}; gmem=${gmem/%G/GB}; gmem=${gmem/%M/MB}

step_route() {   # step -> keywords placed after method/basis
    case $1 in
        opt)       echo "Opt=(MaxCycles=200) Freq SCF=XQC$optscrf" ;;
        optraman)  echo "Opt=(MaxCycles=200) Freq=Raman SCF=XQC$optscrf" ;;
        sp)        echo "SCF=XQC$scrf" ;;
        pop)       echo "Pop=(Full,NPA,Hirshfeld) IOp(6/7=3)$scrf" ;;
        chelpg)    echo "Pop=CHelpG$scrf" ;;
        td)        echo "TD=(NStates=$nstates)$scrf" ;;
        nmr)       if ((spinspin)); then echo "NMR=(GIAO,SpinSpin)$scrf"; else echo "NMR=GIAO$scrf"; fi ;;
        polar)     echo "Polar$scrf" ;;
        epr)       echo "Prop=EPR$scrf" ;;
        wfx)       echo "Output=WFX Density=Current$scrf" ;;
        *) die "unknown step '$1' (see -h)" ;;
    esac
}
for s in ${steps//,/ }; do step_route "$s" >/dev/null; done      # validate early
first=${steps%%,*}
[[ "$first" == opt || "$first" == optraman || "$first" == sp ]] ||
    warn "first step is '$first': it will run on the input geometry, not an optimised one"

# basis set section for Gen: light basis, heavy ECP basis, then the ECP itself
gen_block() {
    echo "$LIGHT_ELEMENTS 0"; echo "$basis"; echo "****"
    echo "$HEAVY_ELEMENTS 0"; echo "$ecp";   echo "****"; echo
    echo "$HEAVY_ELEMENTS 0"; echo "$ecp";   echo
}

# write one --Link1-- step:  step_input NAME STEP CHARGE MULT GEOM(0/1) [chk options]
step_input() {
    local name=$1 step=$2 q=$3 m=$4 withgeom=$5 oldchk=${6:-} chk=${7:-$1.chk} title=${8:-$2}
    local bas="$basis" extra="" gen=0
    # mixed basis: written once in the first step, later steps read it (and
    # the ECP) back from the checkpoint with ChkBasis
    if [[ -n "$HEAVY_ELEMENTS" ]]; then
        if ((withgeom)); then bas="Gen"; extra=" Pseudo=Read"; gen=1; else bas="ChkBasis"; fi
    fi
    [[ -n "$oldchk" ]] && echo "%OldChk=$oldchk"
    echo "%Chk=$chk"
    echo "%NProcShared=$nproc"
    echo "%Mem=$gmem"
    local geo=""; ((withgeom)) || geo=" Geom=Check Guess=Read"
    echo "#P $method/$bas$extra$dispkw $(step_route "$step")$geo"
    echo; echo "$name $title"; echo
    echo "$q $m"
    if ((withgeom)); then cat "$tmp/coords"; fi
    echo
    ((gen)) && gen_block
    [[ "$step" == wfx ]] && { echo "$name.wfx"; echo; }
    return 0
}

tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
mkdir -p "$outdir"
summary="$outdir/prepared.txt"; : > "$summary"
for f in "${files[@]}"; do
    name=$(strip_ext "$f")
    q=$(get_charge "$name" "$cmfile"); m=$(get_mult "$name" "$cmfile")
    check_multiplicity "$f" "$q" "$m"
    metals_first "$f" > "$tmp/coords"
    split_elements "$f"
    [[ -z "$LIGHT_ELEMENTS" ]] && die "$name: Gen basis needs at least one light element"
    [[ "$m" != 1 && ",$steps," == *",nmr,"* ]] && warn "$name: NMR on an open shell molecule (mult $m) is rarely meaningful"
    [[ "$m" == 1 && ",$steps," == *",epr,"* ]] && warn "$name: EPR step on a closed shell molecule gives no hyperfine data"
    d="$outdir/$name"; mkdir -p "$d"
    {
        i=0
        for s in ${steps//,/ }; do
            ((i)) && { echo "--Link1--"; }
            step_input "$name" "$s" "$q" "$m" $((i == 0))
            i=$((i + 1))
        done
    } > "$d/$name.com"
    cmd=("g16 < $name.com > $name.log" "formchk $name.chk $name.fchk")
    if ((fukui)); then
        # N (reference), N+1 (anion) and N-1 (cation) at the same geometry
        if [[ $m == 1 ]]; then mion=2; else mion=$((m - 1)); fi
        {
            step_input "$name" pop "$q" "$m" 0 "$name.chk" "$name-fukui-N.chk" "fukui N (reference)"
            echo "--Link1--"
            step_input "$name" pop "$((q - 1))" "$mion" 0 "$name.chk" "$name-fukui-Nplus1.chk" "fukui N+1 (anion)"
            echo "--Link1--"
            step_input "$name" pop "$((q + 1))" "$mion" 0 "$name.chk" "$name-fukui-Nminus1.chk" "fukui N-1 (cation)"
        } | sed 's/ Guess=Read//' > "$d/$name-fukui.com"      # ions need a fresh guess
        cmd+=("g16 < $name-fukui.com > $name-fukui.log" "for s in N Nplus1 Nminus1; do formchk $name-fukui-\$s.chk $name-fukui-\$s.fchk; done")
    fi
    write_jobscript "g-$name" "$PWD/$d" GAUSSIAN_SETUP "${cmd[@]}" > "$d/run.sh"
    chmod +x "$d/run.sh"
    printf '%-20s q=%-3s mult=%-2s atoms=%-4s heavy=[%s] steps=%s\n' "$name" "$q" "$m" \
        "$(wc -l < "$tmp/coords")" "$HEAVY_ELEMENTS" "$steps" | tee -a "$summary" >&2
done
info "inputs in $outdir/ (summary: $summary). Submit with: bin/qdm-submit $outdir/*/"
