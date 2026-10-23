#!/usr/bin/env bash
# =============================================================================
#  gaussian-ts-input.sh - transition state, IRC and counterpoise jobs
#
#  For every NAME in the list the script looks for these files in -d DIR
#  (any format Open Babel reads):
#     NAME-r.*   reactant (complex)      NAME-ts.*  transition state guess
#     NAME-p.*   product  (complex)
#
#  usage: gaussian-ts-input.sh -t TYPE [options] -l list [-d dir]
#    -t ts      Opt=(TS,CalcFC) from NAME-ts, followed by Freq and an IRC
#       qst3    Opt=(QST3) from NAME-r, NAME-p and NAME-ts (atom order must match)
#       cp      counterpoise corrected interaction energy of NAME-r (or NAME)
#               split in two fragments: -f "1-12" gives the atoms of fragment 1
#    -m METHOD (default B3LYP)  -D DISP (GD3BJ)  -b BASIS (6-311+G(d,p))  -e ECP (LanL2DZ)
#    -s MODEL:SOLVENT           solvent for all steps (default gas phase)
#    -i N       IRC points in each direction (default 30, 0 = no IRC)
#    -f RANGES  fragment 1 atoms for cp, e.g. "1-12,15" (or per molecule in the
#               second column of the list: "name 1-12")
#    -q "Q M Q1 M1 Q2 M2"  charges/multiplicities of complex and fragments for
#               cp (default from chargemult.txt for the complex, fragments 0 1)
#    -c FILE    charge/multiplicity table (default chargemult.txt)
#    -n NPROC  -M MEM  -o OUTDIR (default gaussian-ts)
#
#  Tips kept from the original workflow:
#   * a QST3 run on two or three members of a series gives a TS pattern that
#     can be copied to the others and refined with -t ts
#   * build the guess in a molecular editor by stretching the bonds that break
#     and form; ModRedundant constraints can be added by hand to the route
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"
source "$QDM_ROOT/lib/jobscript.sh"
load_cluster_conf
type=""; method=B3LYP; disp=GD3BJ; basis="6-311+G(d,p)"; ecp=LanL2DZ; solv=""; nirc=30
frag=""; cpq=""; cmfile=chargemult.txt; nproc=${CORES:-8}; mem=${MEM:-16G}; outdir=gaussian-ts
list=""; indir=.
while getopts "t:m:D:b:e:s:i:f:q:c:n:M:o:l:d:h" o; do
    case $o in
        t) type=$OPTARG ;; m) method=$OPTARG ;; D) disp=$OPTARG ;; b) basis=$OPTARG ;;
        e) ecp=$OPTARG ;; s) solv=$OPTARG ;; i) nirc=$OPTARG ;; f) frag=$OPTARG ;;
        q) cpq=$OPTARG ;; c) cmfile=$OPTARG ;; n) nproc=$OPTARG ;; M) mem=$OPTARG ;;
        o) outdir=$OPTARG ;; l) list=$OPTARG ;; d) indir=$OPTARG ;; *) sed -n '2,31p' "$0"; exit 0 ;;
    esac
done
CORES=$nproc; MEM=$mem; JOB_THREADED=1
[[ "$type" =~ ^(ts|qst3|cp)$ ]] || die "choose -t ts, qst3 or cp"
[[ -n "$list" ]] || die "give a list of reaction names with -l"
dispkw=""; [[ "${disp,,}" != none ]] && dispkw=" EmpiricalDispersion=$disp"
scrf=""; [[ -n "$solv" ]] && scrf=" SCRF=(${solv%%:*},Solvent=${solv#*:})"
gmem=${mem%B}; gmem=${gmem/%G/GB}; gmem=${gmem/%M/MB}

# expand "1-3,7" -> "1 2 3 7"
expand_ranges() { echo "$1" | tr ',' '\n' | awk -F- '{ if (NF==2) for(i=$1;i<=$2;i++) print i; else print $1 }' | tr '\n' ' '; }

header() {   # header CHK ROUTE TITLE
    printf '%%Chk=%s\n%%NProcShared=%s\n%%Mem=%s\n#P %s\n\n%s\n\n' "$1" "$nproc" "$gmem" "$2" "$3"
}
basis_kw() { if [[ -n "$HEAVY_ELEMENTS" ]]; then echo "Gen Pseudo=Read"; else echo "$basis"; fi; }
gen_block() {
    [[ -z "$HEAVY_ELEMENTS" ]] && return
    printf '%s 0\n%s\n****\n%s 0\n%s\n****\n\n%s 0\n%s\n\n' "$LIGHT_ELEMENTS" "$basis" \
        "$HEAVY_ELEMENTS" "$ecp" "$HEAVY_ELEMENTS" "$ecp"
}

mkdir -p "$outdir"
while read -r name rest; do
    d="$outdir/$name"; mkdir -p "$d"
    q=$(get_charge "$name" "$cmfile"); m=$(get_mult "$name" "$cmfile")
    case $type in
    ts)
        ts=$(find_structure "$name-ts" "$indir") || { warn "$name: no $name-ts file"; continue; }
        split_elements "$ts"; bk=$(basis_kw)
        {
            header "$name-ts.chk" "$method/$bk$dispkw Opt=(TS,CalcFC,NoEigenTest,MaxCycles=200) Freq SCF=XQC$scrf" "$name TS search"
            echo "$q $m"; metals_first "$ts"; echo; gen_block
            if ((nirc > 0)); then
                echo "--Link1--"
                bk2=$basis; [[ -n "$HEAVY_ELEMENTS" ]] && bk2=ChkBasis
                header "$name-ts.chk" "$method/$bk2$dispkw IRC=(RCFC,MaxPoints=$nirc,StepSize=10) Geom=Check Guess=Read$scrf" "$name IRC from the TS"
                echo "$q $m"; echo
            fi
        } > "$d/$name-ts.com"
        job=$name-ts ;;
    qst3)
        r=$(find_structure "$name-r" "$indir") && p=$(find_structure "$name-p" "$indir") &&
            ts=$(find_structure "$name-ts" "$indir") || { warn "$name: need $name-r, $name-p and $name-ts"; continue; }
        er=$(xyz_body "$r" | awk '{printf $1}'); ep=$(xyz_body "$p" | awk '{printf $1}'); et=$(xyz_body "$ts" | awk '{printf $1}')
        [[ "$er" == "$ep" && "$er" == "$et" ]] || { warn "$name: reactant, product and TS atoms are not in the same order, skipped"; continue; }
        split_elements "$ts"; bk=$(basis_kw)
        {
            header "$name-qst3.chk" "$method/$bk$dispkw Opt=(QST3,CalcFC,MaxCycles=200) Freq SCF=XQC$scrf" "$name reactant"
            echo "$q $m"; metals_first "$r"; echo
            echo "$name product"; echo; echo "$q $m"; metals_first "$p"; echo
            echo "$name TS guess"; echo; echo "$q $m"; metals_first "$ts"; echo
            gen_block
        } > "$d/$name-qst3.com"
        job=$name-qst3 ;;
    cp)
        s=$(find_structure "$name-r" "$indir") || s=$(find_structure "$name" "$indir") || { warn "$name: no structure"; continue; }
        f1=${frag:-$(echo "$rest" | awk '{print $1}')}
        [[ -n "$f1" ]] || { warn "$name: no fragment 1 atoms (-f or third column of the list)"; continue; }
        split_elements "$s"; bk=$(basis_kw)
        qline=${cpq:-"$q $m 0 1 0 1"}
        {
            header "$name-cp.chk" "$method/$bk$dispkw Counterpoise=2 SCF=XQC$scrf" "$name counterpoise, fragment 1 = $f1"
            echo "$qline"
            # original atom numbering is kept here, so the ranges refer to the input file
            xyz_body "$s" | awk -v f1="$(expand_ranges "$f1")" '
                BEGIN{n=split(f1,a," "); for(i=1;i<=n;i++) F[a[i]]=1}
                {printf "%-16s %14.8f %14.8f %14.8f\n", $1 "(Fragment=" ((NR in F)?1:2) ")", $2,$3,$4}'
            echo; gen_block
        } > "$d/$name-cp.com"
        job=$name-cp ;;
    esac
    write_jobscript "g-$job" "$PWD/$d" GAUSSIAN_SETUP "g16 < $job.com > $job.log" "formchk $job.chk $job.fchk" > "$d/run-$type.sh"
    chmod +x "$d/run-$type.sh"
    info "$name: $d/$job.com (submit: qdm-submit -s run-$type.sh $d)"
done < <(sed -e 's/#.*//' -e '/^[[:space:]]*$/d' "$list")
