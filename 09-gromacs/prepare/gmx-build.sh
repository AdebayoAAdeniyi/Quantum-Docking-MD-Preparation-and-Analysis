#!/usr/bin/env bash
# =============================================================================
#  gmx-build.sh - GROMACS topology and solvated box
#
#  Route A (protein-ligand, recommended):  -a AMBER_SYSTEM_DIR
#     converts the tleap system of module 08 (ff14SB + GAFF2 + TIP3P + ions)
#     with ParmEd, so the protein and the ligand come from one consistent force
#     field family, and the same system can be run in Amber and GROMACS.
#  Route B (protein only):  -p protein.pdb [-f amber99sb-ildn] [-w tip3p]
#     pdb2gmx -> editconf (dodecahedron, -b nm) -> solvate -> genion (0.15 M NaCl)
#
#  usage: gmx-build.sh -n NAME (-a DIR | -p PDB) [-f FF] [-w WATER] [-b 1.0] [-s 0.15] [-o OUT]
#  Writes OUT/topol.top, OUT/system.gro, OUT/index.ndx with groups Protein, LIG
#  (if present), Protein_LIG and Water_and_ions, plus position restraint files.
#
#  Why not the original route: the old scripts built the protein with CHARMM36
#  and the ligand with GAFF2 and then edited the ligand 1-4 scaling to make the
#  files combine. Mixing force fields that were parameterised with different
#  rules (charges, 1-4 scaling, water) is not validated; route A avoids it.
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"
name=""; amb=""; pdb=""; ff=amber99sb-ildn; water=tip3p; buf=1.0; salt=0.15; out=""
while getopts "n:a:p:f:w:b:s:o:h" o; do
    case $o in n) name=$OPTARG ;; a) amb=$OPTARG ;; p) pdb=$OPTARG ;; f) ff=$OPTARG ;; w) water=$OPTARG ;; b) buf=$OPTARG ;;
               s) salt=$OPTARG ;; o) out=$OPTARG ;; *) sed -n '2,22p' "$0"; exit 0 ;; esac
done
[[ -n "$name" ]] || die "give -n NAME"
GMX=${GMX:-gmx}; require "${GMX%% *}"
out=${out:-gromacs-systems/$name}; mkdir -p "$out"; out=$(cd "$out" && pwd)
if [[ -n "$amb" ]]; then
    [[ -f "$amb/solv.prmtop" ]] || die "$amb/solv.prmtop not found (run 08-amber/prepare/amber-build.sh)"
    command -v parmed >/dev/null || python3 -c "import parmed" 2>/dev/null || die "ParmEd is needed (it ships with AmberTools)"
    PY=$(command -v parmed >/dev/null && dirname "$(command -v parmed)")/python; [[ -x "$PY" ]] || PY=python3
    "$PY" - "$amb/solv.prmtop" "$amb/solv.inpcrd" "$out" <<'PYEOF'
import sys, parmed
top, crd, out = sys.argv[1:4]
s = parmed.load_file(top, crd)
s.save(f"{out}/topol.top", format="gromacs", overwrite=True)
s.save(f"{out}/system.gro", overwrite=True)
PYEOF
    [[ -s "$out/topol.top" ]] || die "ParmEd conversion failed"
else
    [[ -f "$pdb" ]] || die "give -a AMBER_DIR or -p protein.pdb"
    ( cd "$out" && $GMX pdb2gmx -f "$(realpath "$pdb")" -o protein.gro -p topol.top -ff "$ff" -water "$water" -ignh > pdb2gmx.log 2>&1 ) || die "pdb2gmx failed ($out/pdb2gmx.log)"
    ( cd "$out" && $GMX editconf -f protein.gro -o box.gro -c -d "$buf" -bt dodecahedron > editconf.log 2>&1 &&
      $GMX solvate -cp box.gro -cs spc216.gro -o solv.gro -p topol.top > solvate.log 2>&1 ) || die "box/solvation failed"
    printf 'integrator=steep\nnsteps=0\ncutoff-scheme=Verlet\ncoulombtype=PME\nrcoulomb=1.0\nrvdw=1.0\n' > "$out/ions.mdp"
    ( cd "$out" && $GMX grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 1 > grompp-ions.log 2>&1 &&
      echo SOL | $GMX genion -s ions.tpr -o system.gro -p topol.top -pname NA -nname CL -neutral -conc "$salt" > genion.log 2>&1 ) || die "genion failed"
fi
# index groups: protein, ligand (any non-standard residue that is not water or ion), complex, solvent
lig=$(awk 'NR>2 { r = substr($0, 6, 5); gsub(/ /, "", r); print r }' "$out/system.gro" | sort -u |
      grep -vxE 'ALA|ARG|ASN|ASP|ASH|CYS|CYX|CYM|GLN|GLU|GLH|GLY|HIS|HIE|HID|HIP|ILE|LEU|LYS|LYN|MET|PHE|PRO|SER|THR|TRP|TYR|VAL|ACE|NME|NHE|NALA|NMET|CALA|CSER|SOL|WAT|HOH|NA|CL|Na\+|Cl-|K\+|K|Na|Cl|NA\+|CL-' |
      grep -vE '^[NC](ALA|ARG|ASN|ASP|CYS|GLN|GLU|GLY|HI[SEDP]|ILE|LEU|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL|CYX)$' | grep -v '^[0-9]*$' | head -1)
# default groups first; add only what is missing: LIG, Protein_LIG and Water_and_ions
echo q | $GMX make_ndx -f "$out/system.gro" -o "$out/index.ndx" > "$out/make_ndx.log" 2>&1 || die "make_ndx failed"
grp()  { awk -v g="$1" '$1=="[" && $2==g {found=1} END{exit !found}' "$out/index.ndx"; }
gnum() { awk '$1=="["{n++} END{print n-1}' "$out/index.ndx"; }     # highest existing group number
{
    if [[ -n "$lig" ]] && ! grp LIG; then echo "r $lig"; echo "name $(($(gnum) + 1)) LIG"; fi
    echo q
} > "$out/.ndx.in"
$GMX make_ndx -f "$out/system.gro" -n "$out/index.ndx" -o "$out/index.ndx" < "$out/.ndx.in" >> "$out/make_ndx.log" 2>&1
solute=Protein; grp Protein || solute=non-Water
if [[ -n "$lig" ]] && grp Protein && ! grp Protein_LIG; then
    p=$(awk -v g=Protein '$1=="["{n++} $1=="[" && $2==g {print n-1; exit}' "$out/index.ndx")
    l=$(awk -v g=LIG '$1=="["{n++} $1=="[" && $2==g {print n-1; exit}' "$out/index.ndx")
    { echo "$p | $l"; echo "name $(($(gnum) + 1)) Protein_LIG"; echo q; } > "$out/.ndx.in"
    $GMX make_ndx -f "$out/system.gro" -n "$out/index.ndx" -o "$out/index.ndx" < "$out/.ndx.in" >> "$out/make_ndx.log" 2>&1
    solute=Protein_LIG
fi
if ! grp Water_and_ions; then
    sg=$(awk -v g="$solute" '$1=="["{n++} $1=="[" && $2==g {print n-1; exit}' "$out/index.ndx")
    { echo "! $sg"; echo "name $(($(gnum) + 1)) Water_and_ions"; echo q; } > "$out/.ndx.in"
    $GMX make_ndx -f "$out/system.gro" -n "$out/index.ndx" -o "$out/index.ndx" < "$out/.ndx.in" >> "$out/make_ndx.log" 2>&1
fi
rm -f "$out/.ndx.in"
natoms=$(sed -n 2p "$out/system.gro" | tr -d ' ')
info "$name: $natoms atoms, ligand residue: ${lig:-none}; files in $out (topol.top, system.gro, index.ndx)"
