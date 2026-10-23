#!/usr/bin/env bash
# =============================================================================
#  amber-ligand-params.sh - GAFF2 parameters with AM1-BCC charges for ligands
#
#  For each ligand (any format Open Babel reads, hydrogens included):
#    antechamber  -> NAME/LIG.mol2  (GAFF2 atom types, AM1-BCC charges)
#    parmchk2     -> NAME/LIG.frcmod (missing parameters)
#  The net charge comes from chargemult.txt (or the formal charges in the file).
#  Every frcmod line that parmchk2 had to guess ("ATTN, need revision") is
#  reported, because those are the parameters to check by hand or refit.
#
#  usage: amber-ligand-params.sh [-c chargemult.txt] [-r RES] [-m bcc|gas|resp-file]
#                                [-o outdir] [-l list -d dir | files...]
#    -r  three-letter residue name (default: LIG, or the first 3 characters of
#        the name in upper case when several ligands are given)
#    -m  charge model: bcc (AM1-BCC, default) or gas (Gasteiger, fast, for tests)
#  RESP charges from a Gaussian HF/6-31G* ESP job (the R.E.D. route of the
#  original scripts) can be used by giving the Gaussian log as input with -m resp.
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"
cmfile=chargemult.txt; res=""; model=bcc; out=amber-ligands; list=""; indir=.
while getopts "c:r:m:o:l:d:h" o; do
    case $o in c) cmfile=$OPTARG ;; r) res=$OPTARG ;; m) model=$OPTARG ;; o) out=$OPTARG ;; l) list=$OPTARG ;; d) indir=$OPTARG ;;
               *) sed -n '2,19p' "$0"; exit 0 ;; esac
done
shift $((OPTIND-1))
require antechamber parmchk2
mapfile -t files < <(gather_structures "$list" "$indir" "$@")
((${#files[@]})) || die "no ligands"
for f in "${files[@]}"; do
    n=$(strip_ext "$f"); d="$out/$n"; mkdir -p "$d"
    r=${res:-$( ((${#files[@]} > 1)) && echo "${n:0:3}" | tr a-z A-Z || echo LIG )}
    q=$(get_charge "$n" "$cmfile")
    if [[ ! -f "$cmfile" ]] || ! grep -qw "^$n" "$cmfile"; then     # fall back to formal charges
        fq=$(obabel "$f" -osmi 2>/dev/null | head -1 | awk '{s=$1; n=gsub(/\+\]/,"",s); m=gsub(/-\]/,"",s); print n-m}')
        [[ -n "$fq" && "$fq" != 0 ]] && { q=$fq; info "$n: net charge $q from the formal charges"; }
    fi
    to_format "$f" "$d/input.mol2"
    nH=$(awk '/@<TRIPOS>ATOM/{a=1;next} /@<TRIPOS>/{a=0} a && $6 ~ /^H/' "$d/input.mol2" | wc -l)
    ((nH > 0)) || warn "$n: no hydrogens in the structure; add them for the right pH first (convert-structures.sh -p 7.4)"
    case $model in
        bcc|gas) (cd "$d" && antechamber -i input.mol2 -fi mol2 -o LIG.mol2 -fo mol2 -c "$model" -nc "$q" -at gaff2 -rn "$r" -pf y > antechamber.log 2>&1) ;;
        resp)    (cd "$d" && antechamber -i "$(realpath "$f")" -fi gout -o LIG.mol2 -fo mol2 -c resp -nc "$q" -at gaff2 -rn "$r" -pf y > antechamber.log 2>&1) ;;
        *) die "unknown charge model $model" ;;
    esac
    [[ -s "$d/LIG.mol2" ]] || { warn "$n: antechamber failed, see $d/antechamber.log"; continue; }
    (cd "$d" && parmchk2 -i LIG.mol2 -f mol2 -o LIG.frcmod -s gaff2)
    qsum=$(awk '/@<TRIPOS>ATOM/{a=1;next} /@<TRIPOS>/{a=0} a{s+=$9} END{printf "%.3f", s}' "$d/LIG.mol2")
    nattn=$(grep -c "ATTN" "$d/LIG.frcmod" || true)
    info "$n: residue $r, net charge $q (sum of partial charges $qsum), $model charges, $nattn guessed parameters"
    ((nattn)) && grep "ATTN" "$d/LIG.frcmod" | sed 's/^/      /' >&2
done
