#!/usr/bin/env bash
# =============================================================================
#  gamess-input.sh - GAMESS-US and Firefly inputs, including NBO natural energy
#                    decomposition (NEDA) between two or more fragments
#
#  The basis set of every atom is written explicitly in $DATA and the ECPs in
#  $ECP, both taken from the Basis Set Exchange library (pip install
#  basis_set_exchange), so light atoms can use one basis and metals another
#  (the original scripts kept hand-made library files for this).
#
#  usage: gamess-input.sh [options] [-l list -d dir | files...]
#    -P PROG     gamess (GAMESS-US, default) or firefly
#    -j JOB      energy | opt | hessian | optfreq | neda      (default energy)
#    -m METHOD   DFT functional (B3LYP, PBE0, M06 ...) or HF  (default B3LYP)
#    -b BASIS    basis for light atoms, any BSE name     (default 6-31G*)
#    -e ECP      basis + ECP for heavy atoms              (default LANL2DZ)
#    -f GROUPS   NEDA fragments as atom ranges separated by "/",
#                e.g. "1-12/13-30" (or per molecule in the list: "name 1-12/13-30");
#                without -f, NAME-f1.*, NAME-f2.* ... fragment files are used to
#                find which atoms of NAME belong to each fragment
#    -c FILE     charge/multiplicity table (default chargemult.txt)
#    -n NPROC -M MEM -o OUTDIR (default gamess-jobs)
#  NEDA needs GAMESS or Firefly linked to NBO 5/6/7. The keyword "NAO NEDA"
#  with atom lists is used, as in the NBO NEDA tutorial
#  (https://nbo6.chem.wisc.edu/tut_neda.htm).
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"
source "$QDM_ROOT/lib/jobscript.sh"
load_cluster_conf
prog=gamess; job=energy; method=B3LYP; basis="6-31G*"; ecp=LANL2DZ; frags=""; cmfile=chargemult.txt
nproc=${CORES:-8}; mem=${MEM:-16G}; outdir=gamess-jobs; list=""; indir=.
while getopts "P:j:m:b:e:f:c:n:M:o:l:d:h" o; do
    case $o in
        P) prog=${OPTARG,,} ;; j) job=$OPTARG ;; m) method=${OPTARG^^} ;; b) basis=$OPTARG ;; e) ecp=$OPTARG ;;
        f) frags=$OPTARG ;; c) cmfile=$OPTARG ;; n) nproc=$OPTARG ;; M) mem=$OPTARG ;; o) outdir=$OPTARG ;;
        l) list=$OPTARG ;; d) indir=$OPTARG ;; *) sed -n '2,26p' "$0"; exit 0 ;;
    esac
done
shift $((OPTIND-1))
[[ "$prog" =~ ^(gamess|firefly)$ ]] || die "-P must be gamess or firefly"
require bse
CORES=$nproc; MEM=$mem
memmb=$(echo "$mem" | awk '{v=$1+0; u=toupper($1); if (u ~ /T/) v*=1048576; else if (u ~ /G/) v*=1024; print v}')
mwords=$(( memmb * 3 / 4 / 8 / nproc ))          # 1 word = 8 bytes, per process, 25 % headroom
cache="${QDM_CACHE:-$HOME/.cache/qdm}/bse"; mkdir -p "$cache"

# ---- basis set helpers (BSE, GAMESS-US format) ----
bse_raw() {   # bse_raw BASIS ELEMENT -> cached text
    local f="$cache/$(echo "$1" | tr '/*+()' '_spLR')-$2.gms"
    [[ -s "$f" ]] || bse get-basis "$1" gamess_us --elements "$2" > "$f" 2>/dev/null || { rm -f "$f"; return 1; }
    echo "$f"
}
shells() {    # shells BASIS EL : shell lines of EL (the block after its full name)
    local f; f=$(bse_raw "$1" "$2") || die "basis $1 not available for $2 in the Basis Set Exchange"
    awk '/^\$DATA/{d=1; next} d && /^[A-Z]+$/ {inb=1; next} inb && /^$/ {exit} inb && !/^\$END/ {print "  " $0}' "$f"
}
ecp_block() { # ecp_block BASIS EL : "EL-ECP GEN n l" block, empty when the basis has no ECP for EL
    local f; f=$(bse_raw "$1" "$2") || return 0
    awk '/^\$ECP/{e=1; next} e && /^\$END/ {exit} e' "$f"
}
expand() { echo "$1" | tr ',' '\n' | awk -F- 'NF==2{for(i=$1;i<=$2;i++) printf "%d,", i; next} {printf "%d,", $1}' | sed 's/,$//'; }

mkdir -p "$outdir"
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
while read -r name rest; do
    f=$(find_structure "$name" "$indir") || { warn "no structure for $name"; continue; }
    q=$(get_charge "$name" "$cmfile"); m=$(get_mult "$name" "$cmfile")
    check_multiplicity "$f" "$q" "$m"
    xyz_body "$f" > "$tmp/xyz"; split_elements "$f"
    d="$outdir/$name"; mkdir -p "$d"; inp="$d/$name-$job.inp"
    scf=RHF; ((m > 1)) && scf=UHF
    case $job in
        energy) rt=ENERGY ;; opt) rt=OPTIMIZE ;; hessian) rt=HESSIAN ;; optfreq) rt=OPTIMIZE ;; neda) rt=ENERGY ;;
        *) die "unknown job $job" ;;
    esac
    dft=""; [[ "$method" != HF ]] && dft=" DFTTYP=$method"
    pp=""; if [[ -n "$HEAVY_ELEMENTS" ]]; then [[ $prog == gamess ]] && pp=" PP=READ" || pp=" ECP=READ"; fi
    sph=""; [[ "${basis,,}" =~ ^(def2|cc-|aug-cc) ]] && sph=" ISPHER=1"
    {
        echo " \$CONTRL SCFTYP=$scf RUNTYP=$rt ICHARG=$q MULT=$m$dft$pp$sph"
        echo "         COORD=UNIQUE UNITS=ANGS MAXIT=200 NOSYM=1 \$END"
        if [[ $prog == gamess ]]; then echo " \$SYSTEM MWORDS=$mwords \$END"
        else echo " \$SYSTEM MWORDS=$mwords TIMLIM=100000 \$END"; echo " \$P2P P2P=.T. DLB=.T. \$END"; fi
        echo " \$SCF DIRSCF=.T. DIIS=.T. DAMP=.T. \$END"
        [[ $rt == OPTIMIZE ]] && echo " \$STATPT NSTEP=200 $([[ $job == optfreq ]] && echo 'HSSEND=.T.') \$END"
        if [[ $job == neda ]]; then
            fr=${frags:-$(echo "$rest" | awk '{print $1}')}
            if [[ -z "$fr" ]]; then     # fragment files NAME-f1, NAME-f2 ...: match atoms by position
                fr=""
                for ff in $(ls "$indir/$name"-f[0-9]* 2>/dev/null | sort -V); do
                    idx=$(xyz_body "$ff" | awk 'NR==FNR{x[NR]=$2" "$3" "$4; e[NR]=$1; n=NR; next}
                          { for(i=1;i<=n;i++){ split(x[i],c," "); if ($1==e[i] && (c[1]-$2)^2+(c[2]-$3)^2+(c[3]-$4)^2 < 1e-4) printf "%d,", i } }' "$tmp/xyz" - | sed 's/,$//')
                    fr+="${fr:+/}$idx"
                done
            fi
            [[ "$fr" == */* ]] || die "$name: NEDA needs two or more fragments (-f \"1-12/13-30\" or $name-f1, $name-f2 files)"
            groups=""; for g in ${fr//\// }; do groups+=" ($(expand "$g"))"; done
            ntot=$(wc -l < "$tmp/xyz"); nin=$(echo "$groups" | grep -o '[0-9]\+' | wc -l)
            ((nin == ntot)) || warn "$name: fragments cover $nin of $ntot atoms"
            echo " \$NBO \$END"
            echo " \$DEL"
            if ((m > 1)); then echo " ALPHA NAO NEDA$groups END"; echo " BETA NAO NEDA$groups END"
            else echo " NAO NEDA$groups END"; fi
            echo " \$END"
        fi
        echo " \$DATA"; echo "$name $job"; echo "C1"
        while read -r el x y z; do
            z_=$(echo "$el" | awk -f "$QDM_ROOT/lib/awk/atomic_number.awk")
            printf ' %-3s %5.1f %14.8f %14.8f %14.8f\n' "$el" "$z_" "$x" "$y" "$z"
            if is_heavy "$el"; then shells "$ecp" "$el"; else shells "$basis" "$el"; fi
            echo
        done < "$tmp/xyz"
        echo " \$END"
        if [[ -n "$HEAVY_ELEMENTS" ]]; then
            echo " \$ECP"
            declare -A seen=()
            while read -r el _; do
                u=${el^^}
                if ! is_heavy "$el"; then echo "  $u-ECP NONE"
                elif [[ -n "${seen[$u]:-}" ]]; then echo "  $u-ECP"          # same potential as before
                else seen[$u]=1; ecp_block "$ecp" "$el" | sed 's/^/  /'; fi
            done < "$tmp/xyz"
            unset seen
            echo " \$END"
        fi
    } > "$inp"
    if [[ $prog == gamess ]]; then
        cmd="\${GAMESS_EXE:-rungms} $name-$job.inp 00 $nproc > $name-$job.log 2>&1"; setup=GAMESS_SETUP
    else
        cmd="\${FIREFLY_EXE:-firefly} -r -f -p -stdext -i $name-$job.inp -o $name-$job.log"; setup=FIREFLY_SETUP
    fi
    write_jobscript "${prog:0:2}-$name" "$PWD/$d" "$setup" "$cmd" > "$d/run-$job.sh"; chmod +x "$d/run-$job.sh"
    info "$name ($prog, $job): $inp"
done < <( { [[ -n "$list" ]] && sed -e 's/#.*//' -e '/^[[:space:]]*$/d' "$list"; for a in "$@"; do strip_ext "$a"; done; } )
