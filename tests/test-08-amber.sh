#!/usr/bin/env bash
# Amber: real ligand parameterisation and a real solvated build (small, so it is quick).
source "$(dirname "$0")/testlib.sh"
need antechamber; need tleap
A=$QDM_ROOT/08-amber
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT; cd "$W"
echo "== ligand parameters (GAFF2, Gasteiger charges for speed)"
fetch "https://raw.githubusercontent.com/ccsb-scripps/AutoDock-Vina/develop/example/basic_docking/data/1iep_ligand.sdf" 1iep_ligand.sdf || exit 0
cp "$DATA/1iep_ligand.sdf" STI.sdf
"$A/prepare/amber-ligand-params.sh" -m gas -o lig STI.sdf 2>/dev/null
check "mol2 written"           "$([[ -s lig/STI/LIG.mol2 ]] && echo yes)" yes
check "net charge +1 detected" "$(awk '/@<TRIPOS>ATOM/{a=1;next} /@<TRIPOS>/{a=0} a{s+=$9} END{printf "%.2f", s}' lig/STI/LIG.mol2)" 1.00 0.01
check "GAFF2 atom types"       "$(awk '/@<TRIPOS>ATOM/{a=1;next} /@<TRIPOS>/{a=0} a && $6=="ca"' lig/STI/LIG.mol2 | wc -l | awk '{print ($1>5)?"yes":"no"}')" yes
check "frcmod written"         "$([[ -s lig/STI/LIG.frcmod ]] && echo yes)" yes
echo "== solvated ligand box"
"$A/prepare/amber-build.sh" -n lw -L lig/STI -b 8 -s 0.15 -o sys 2>/dev/null
check "solvated topology"      "$([[ -s sys/solv.prmtop ]] && echo yes)" yes
check "dry topologies for MM-PBSA" "$(ls sys/complex.prmtop sys/ligand.prmtop 2>/dev/null | wc -l)" 2
nw=$(grep -c "O   WAT" sys/solv.pdb)
check "waters added"           "$(( nw > 500 ? 1 : 0 ))" 1
check "salt added"             "$(grep -c "Na+\|Cl-" sys/solv.pdb | awk '{print ($1>0)?"yes":"no"}')" yes
echo "== MD inputs"
"$A/run/amber-md-inputs.sh" -d sys -T 300 -H 10 -D 10 -E 10 -P 0.01 -f 1 2>/dev/null
check "six stages"             "$(ls sys/{min1,min2,heat,density,equil,prod}.in 2>/dev/null | wc -l)" 6
check "temperature ramp"       "$(grep -c "TEMP0" sys/heat.in)" 2
check "MC barostat in NPT"     "$(grep -c "barostat=2" sys/density.in)" 1
check "restart logic"          "$(grep -c 'already done' sys/run-md.sh)" 1
if [[ -n "${QDM_SLOW:-}" ]] && command -v sander >/dev/null; then
    sed -i 's/maxcyc=5000, ncyc=2500/maxcyc=200, ncyc=100/' sys/min1.in sys/min2.in
    (cd sys && sander -O -i min1.in -p solv.prmtop -c solv.inpcrd -o min1.out -r min1.rst7 -ref solv.inpcrd > /dev/null 2>&1)
    check "minimisation runs"  "$([[ -s sys/min1.rst7 ]] && echo yes)" yes
fi
finish
