#!/usr/bin/env bash
# =============================================================================
#  build-3d-from-smiles.sh - SMILES table -> optimised 3D structures
#
#  Input: a text file with "name SMILES" per line, for example
#      1-H    Oc1ccccc1C(=O)c1ccccc1
#      2-OMe  COc1ccc(C(=O)c2ccccc2)c(O)c1
#
#  For every molecule: protonate (optional pH), generate 3D coordinates, search
#  conformers with a genetic algorithm and minimise with MMFF94 (falls back to
#  UFF for elements MMFF94 does not cover). The lowest-energy conformer is
#  written in every requested format.
#
#  This replaces download-2d-nist-optimized-to-3d-obabel, which needed internet
#  access to the NIST WebBook and ran three obminimize passes by hand.
#
#  usage: build-3d-from-smiles.sh -i smiles.txt [-o outdir] [-f xyz,pdb,mol2]
#                                 [-p 7.4] [-n 30] [-s 2000]
#    -p  pH for protonation (default: keep the SMILES as written)
#    -n  number of conformers to try (default 30)
#    -s  minimisation steps (default 2000)
# =============================================================================
source "$(dirname "$0")/../lib/common.sh"
require obabel
outdir=structures; formats=xyz,pdb,mol2; ph=""; nconf=30; steps=2000; input=""
while getopts "i:o:f:p:n:s:h" o; do
    case $o in i) input=$OPTARG ;; o) outdir=$OPTARG ;; f) formats=$OPTARG ;;
               p) ph=$OPTARG ;; n) nconf=$OPTARG ;; s) steps=$OPTARG ;;
               *) sed -n '2,25p' "$0"; exit 0 ;; esac
done
[[ -f "$input" ]] || die "give the SMILES table with -i (see -h)"
mkdir -p "$outdir"
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT

while read -r name smiles _; do
    [[ -z "$name" || "$name" == \#* ]] && continue
    info "$name  $smiles"
    echo "$smiles $name" > "$tmp/in.smi"
    phopt=(); [[ -n "$ph" ]] && phopt=(-p "$ph")
    # 3D build + quick clean-up
    obabel "$tmp/in.smi" -O "$tmp/raw.sdf" --gen3d "${phopt[@]}" 2>/dev/null
    ff=MMFF94
    obminimize -ff MMFF94 -n 1 "$tmp/raw.sdf" >/dev/null 2>&1 || ff=UFF
    # conformer search (genetic algorithm), then a full minimisation
    obabel "$tmp/raw.sdf" -O "$tmp/conf.sdf" --conformer --nconf "$nconf" --score energy --ff "$ff" >/dev/null 2>&1 ||
        cp "$tmp/raw.sdf" "$tmp/conf.sdf"
    obminimize -ff "$ff" -sd -c 1e-6 -n "$steps" -osdf "$tmp/conf.sdf" > "$tmp/min.sdf" 2>/dev/null
    [[ -s "$tmp/min.sdf" ]] || die "minimisation failed for $name"
    energy=$(obenergy -ff "$ff" "$tmp/min.sdf" 2>/dev/null | awk '/TOTAL ENERGY/{print $4, $5}')
    for f in ${formats//,/ }; do
        obabel "$tmp/min.sdf" -O "$outdir/$name.$f" --title "$name" 2>/dev/null
    done
    printf '%-14s %-8s %s\n' "$name" "$ff" "$energy" >> "$outdir/build-energies.txt"
done < "$input"
info "structures written to $outdir/ (force-field energies in build-energies.txt)"
