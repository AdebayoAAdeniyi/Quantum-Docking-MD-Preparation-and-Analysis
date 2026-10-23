#!/usr/bin/env bash
# =============================================================================
#  ligand-pdb-for-md.sh - give a ligand PDB unique atom names (C1 C2 ... O1 H1)
#  and one residue name, as antechamber/tleap/pdb2gmx expect.
#
#  Converts the input first, so xyz, mol2 or sdf work as well as pdb.
#  (Replaces numbering-pdb / numbering-2-pdb.)
#
#  usage: ligand-pdb-for-md.sh INPUT [RESNAME=LIG] > ligand.pdb
# =============================================================================
source "$(dirname "$0")/../lib/common.sh"
in="${1:?usage: ligand-pdb-for-md.sh INPUT [RESNAME]}"; res="${2:-LIG}"
tmp=$(mktemp --suffix=.pdb); trap 'rm -f "$tmp"' EXIT
to_format "$in" "$tmp"
awk -v res="${res:0:3}" '
/^(ATOM|HETATM)/ {
    el = substr($0, 77, 2); gsub(/ /, "", el)
    if (el == "") { el = substr($0, 13, 4); gsub(/[^A-Za-z]/, "", el); el = substr(el, 1, 1) }
    n[el]++; name = el n[el]
    # PDB rule: one-letter elements start in column 14, two-letter in column 13
    fname = (length(el) == 1 && length(name) < 4) ? " " name : name
    printf "HETATM%5d %-4s %3s A   1    %s%6.2f%6.2f          %2s\n", ++i, fname, res, substr($0, 31, 24), 1.0, 0.0, el
    next }
/^CONECT/ { print; next }
END { print "END" }' "$tmp"
