#!/usr/bin/env bash
# =============================================================================
#  amber-build.sh - solvated protein-ligand (or ligand-only, or protein-only)
#                   systems with tleap, plus the dry topologies for MM-PBSA
#
#  usage: amber-build.sh -n NAME [-p receptor.pdb] [-L ligdir ...] [options]
#    -p PDB     receptor; cleaned with pdb4amber (waters and hydrogens removed,
#               histidine states kept, TER cards at chain breaks)
#    -L DIR     ligand folder made by amber-ligand-params.sh (LIG.mol2, LIG.frcmod);
#               repeat -L for several ligands or cofactors
#    -w MODEL   water: tip3p (default), opc, spce
#    -b BUF     distance from solute to box edge in A (default 10)
#    -B oct|box truncated octahedron (default) or rectangular box
#    -s M       salt concentration in mol/L added on top of neutralisation (default 0.15)
#    -X         keep crystallographic waters (HOH within the PDB)
#    -o DIR     output folder (default amber-systems/NAME)
#  Force fields: ff14SB (protein), GAFF2 (ligands), Joung-Cheatham ions for the
#  chosen water. Writes: solv.prmtop/.inpcrd/.pdb, and for MM-PBSA
#  complex.prmtop, receptor.prmtop, ligand.prmtop (dry), plus leap.in and leap.log.
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"
name=""; rec=""; ligs=(); water=tip3p; buf=10; shape=oct; salt=0.15; xtal=0; out=""
while getopts "n:p:L:w:b:B:s:Xo:h" o; do
    case $o in n) name=$OPTARG ;; p) rec=$OPTARG ;; L) ligs+=("$OPTARG") ;; w) water=$OPTARG ;; b) buf=$OPTARG ;;
               B) shape=$OPTARG ;; s) salt=$OPTARG ;; X) xtal=1 ;; o) out=$OPTARG ;; *) sed -n '2,20p' "$0"; exit 0 ;; esac
done
[[ -n "$name" ]] || die "give a system name with -n"
[[ -n "$rec" || ${#ligs[@]} -gt 0 ]] || die "nothing to build (-p receptor and/or -L ligand)"
require tleap
out=${out:-amber-systems/$name}; mkdir -p "$out"
case $water in
    tip3p) wbox=TIP3PBOX; wrc=leaprc.water.tip3p ;;
    opc)   wbox=OPCBOX;   wrc=leaprc.water.opc ;;
    spce)  wbox=SPCBOX;   wrc=leaprc.water.spce ;;
    *) die "water model $water not handled" ;;
esac
if [[ -n "$rec" ]]; then
    require pdb4amber
    opts=(--nohyd); ((xtal)) || opts+=(--dry)
    pdb4amber -i "$rec" -o "$out/receptor-clean.pdb" "${opts[@]}" > "$out/pdb4amber.log" 2>&1 || die "pdb4amber failed, see $out/pdb4amber.log"
    grep -i "gap" "$out/pdb4amber.log" | sed 's/^/  /' >&2 || true
fi
{
    echo "source leaprc.protein.ff14SB"
    echo "source leaprc.gaff2"
    echo "source $wrc"
    echo "set default PBRadii mbondi2"
    parts=()
    i=0
    for L in "${ligs[@]}"; do
        i=$((i + 1)); [[ -f "$L/LIG.mol2" && -f "$L/LIG.frcmod" ]] || die "$L has no LIG.mol2/LIG.frcmod"
        echo "loadAmberParams $(realpath "$L/LIG.frcmod")"
        echo "L$i = loadMol2 $(realpath "$L/LIG.mol2")"
        echo "check L$i"
        parts+=("L$i")
    done
    if ((${#ligs[@]})); then
        if ((${#ligs[@]} > 1)); then echo "LIG = combine { ${parts[*]} }"; else echo "LIG = L1"; fi
    fi
    [[ -n "$rec" ]] && echo "REC = loadPdb $(realpath "$out/receptor-clean.pdb")"
    if [[ -n "$rec" && ${#ligs[@]} -gt 0 ]]; then echo "COM = combine { REC LIG }"
    elif [[ -n "$rec" ]]; then echo "COM = copy REC"
    else echo "COM = copy LIG"; fi
    [[ -n "$rec" ]]       && echo "saveAmberParm REC receptor.prmtop receptor.inpcrd"
    ((${#ligs[@]}))       && echo "saveAmberParm LIG ligand.prmtop ligand.inpcrd"
    echo "saveAmberParm COM complex.prmtop complex.inpcrd"
    echo "savePdb COM complex-dry.pdb"
    echo "charge COM"
    if [[ $shape == oct ]]; then echo "solvateOct COM $wbox $buf"; else echo "solvateBox COM $wbox $buf"; fi
    echo "addIons2 COM Na+ 0"
    echo "addIons2 COM Cl- 0"
    echo "saveAmberParm COM solv.prmtop solv.inpcrd"
    echo "savePdb COM solv.pdb"
    echo "quit"
} > "$out/leap-1.in"
# first pass: neutralise and count waters, so the salt can be added at the requested molarity
(cd "$out" && tleap -f leap-1.in > leap-1.log 2>&1) || die "tleap failed, see $out/leap-1.log"
grep -q "Exiting LEaP: Errors = 0" "$out/leap-1.log" || { grep -i "error\|fatal" "$out/leap-1.log" | head -5 >&2; die "tleap reported errors ($out/leap-1.log)"; }
nwat=$(grep -c "^ATOM.*O   WAT" "$out/solv.pdb" || true)
# ions per pair = molarity x water volume, with 55.5 mol/L water
npair=$(awk -v s="$salt" -v w="$nwat" 'BEGIN{printf "%d", s * w / 55.5 + 0.5}')
sed "s/addIons2 COM Cl- 0/addIons2 COM Cl- 0\naddIonsRand COM Na+ $npair Cl- $npair/" "$out/leap-1.in" > "$out/leap.in"
(cd "$out" && tleap -f leap.in > leap.log 2>&1) || die "tleap failed, see $out/leap.log"
grep -q "Exiting LEaP: Errors = 0" "$out/leap.log" || die "tleap reported errors ($out/leap.log)"
rm -f "$out/leap-1.in" "$out/leap-1.log" "$out/leap.log.1"
natom=$(awk '/%FLAG POINTERS/{getline; getline; print $1; exit}' "$out/solv.prmtop")
qtot=$(grep -m1 "Total unperturbed charge" "$out/leap.log" | awk '{print $NF}')
warns=$(grep -o "Warnings = [0-9]*" "$out/leap.log" | tail -1 | awk "{print \$3}")
info "$name: $natom atoms, $nwat waters, +$npair Na+/Cl- pairs (~$salt M), solute charge $qtot, $warns tleap warnings (see $out/leap.log)"
