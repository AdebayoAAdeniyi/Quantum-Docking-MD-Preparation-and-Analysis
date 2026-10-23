#!/usr/bin/env bash
# =============================================================================
#  split-multistructure.sh - split a multi-molecule file (docking poses, a
#  ZINC/PubChem download, an NMR ensemble, a trajectory PDB) into one file per
#  molecule, named after the molecule title (e.g. ZINC08445875) or MODEL number.
#  Repeated titles get _01, _02 ... so no pose is overwritten.
#
#  Replaces splitsdf.csh, splitmol2.csh, splitpdb.csh and the MOE/Glide
#  splitters, which each did the same job with csplit for one format.
#
#  usage: split-multistructure.sh FILE [-o outdir] [-f outformat] [-k N]
#    -f  output format (default: same as input)
#    -k  keep only the first N structures per title (e.g. -k 1 = best pose)
# =============================================================================
source "$(dirname "$0")/../lib/common.sh"
require obabel
in="${1:-}"; [[ -f "$in" ]] || { sed -n '2,15p' "$0"; exit 1; }; shift
outdir=split; ofmt=""; keep=0
while getopts "o:f:k:" o; do
    case $o in o) outdir=$OPTARG ;; f) ofmt=$OPTARG ;; k) keep=$OPTARG ;; *) exit 1 ;; esac
done
ifmt=$(ext_of "$in"); ofmt=${ofmt:-$ifmt}
mkdir -p "$outdir"
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
obabel -i"$ifmt" "$in" -o"$ifmt" -O "$tmp/part.$ifmt" -m 2>/dev/null
declare -A count
n=0
for p in $(ls "$tmp"/part*."$ifmt" | sort -V); do
    raw=$(obabel -i"$ifmt" "$p" -otxt 2>/dev/null | head -n 1)
    # prefer a database identifier hidden in long titles ("ZINC123.mol2[pose 3]")
    title=$(grep -oE 'ZINC[0-9]+|CHEMBL[0-9]+|CID[0-9]+' <<< "$raw" | head -n 1 || true)
    [[ -z "$title" ]] && title=$(tr -cs 'A-Za-z0-9_.-' '_' <<< "$raw" | sed 's/^_*//; s/_*$//')
    [[ -z "$title" ]] && title="mol$((n+1))"
    c=$(( ${count[$title]:-0} + 1 )); count[$title]=$c
    (( keep > 0 && c > keep )) && continue
    out="$outdir/$title"; (( c > 1 )) && out=$(printf '%s_%02d' "$out" $((c-1)))
    obabel -i"$ifmt" "$p" -o"$ofmt" -O "$out.$ofmt" 2>/dev/null
    n=$((n+1))
done
info "$n structures written to $outdir/"
