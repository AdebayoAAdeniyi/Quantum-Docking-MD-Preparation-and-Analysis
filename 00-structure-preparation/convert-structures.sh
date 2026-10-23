#!/usr/bin/env bash
# =============================================================================
#  convert-structures.sh - convert any set of structures to the formats needed
#                          by the QM, MD and docking workflows
#
#  Accepts every format Open Babel reads (xyz, pdb, mol2, sdf, mol, cml, smi,
#  Gaussian .log/.com, ...). SMILES and flat 2D files get 3D coordinates.
#
#  usage: convert-structures.sh [-l list | files...] [-d indir] [-o outdir]
#                               [-f xyz,pdb,mol2] [-H] [-p pH] [-c]
#    -l  list of molecule names (found in -d as name.xyz/.pdb/.mol2/...)
#    -H  add hydrogens          -p  protonate for this pH (implies -H)
#    -c  centre the molecule at the origin
#  example: convert-structures.sh -l list -d raw -o prepared -f xyz,pdb -H
# =============================================================================
source "$(dirname "$0")/../lib/common.sh"
require obabel
list=""; indir=.; outdir=converted; formats=xyz,pdb,mol2; extra=()
while getopts "l:d:o:f:Hp:ch" o; do
    case $o in l) list=$OPTARG ;; d) indir=$OPTARG ;; o) outdir=$OPTARG ;; f) formats=$OPTARG ;;
               H) extra+=(-h) ;; p) extra+=(-p "$OPTARG") ;; c) extra+=(-c) ;;
               *) sed -n '2,17p' "$0"; exit 0 ;; esac
done
shift $((OPTIND-1))
files=("$@")
if [[ -n "$list" ]]; then
    while read -r n; do
        f=$(find_structure "$n" "$indir") || { warn "no structure for $n in $indir"; continue; }
        files+=("$f")
    done < <(read_list "$list")
fi
((${#files[@]})) || die "nothing to convert (give -l list or file names)"
mkdir -p "$outdir"
for f in "${files[@]}"; do
    name=$(strip_ext "$f")
    for fmt in ${formats//,/ }; do
        to_format "$f" "$outdir/$name.$fmt" "${extra[@]}"
    done
    info "$name: $(natoms "$outdir/$name.${formats%%,*}") atoms -> ${formats}"
done
