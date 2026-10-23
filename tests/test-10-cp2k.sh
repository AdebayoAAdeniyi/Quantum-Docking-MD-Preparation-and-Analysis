#!/usr/bin/env bash
# CP2K: a real short MM run on an Amber system, and the QM/MM input checked.
source "$(dirname "$0")/testlib.sh"
need tleap; need antechamber; need cpptraj
C=$QDM_ROOT/10-cp2k; A=$QDM_ROOT/08-amber
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT; cd "$W"
# a small solvated system: methanol in water
printf 'CO methanol\n' > m.smi
"$A/prepare/amber-ligand-params.sh" -m gas -o lig m.smi 2>/dev/null
check "ligand parameters"   "$([[ -s lig/m/LIG.mol2 ]] && echo yes)" yes
"$A/prepare/amber-build.sh" -n mw -L lig/m -b 7 -B box -o sys 2>/dev/null
check "solvated system"     "$([[ -s sys/solv.prmtop ]] && echo yes)" yes
cp sys/solv.inpcrd sys/prod.rst7 2>/dev/null
"$C/cp2k-input.sh" -d sys -j mm -n 5 -E NVT 2>/dev/null
check "MM input written"    "$([[ -s sys/sys-mm.inp ]] && echo yes)" yes
check "Amber 1-4 scaling"   "$(grep -c 'EI_SCALE14 0.8333333' sys/sys-mm.inp)" 1
check "elements declared"   "$(grep -c '&KIND' sys/sys-mm.inp | awk '{print ($1>5)?"yes":"no"}')" yes
check "PDB coordinates"     "$(grep -c 'COORD_FILE_FORMAT PDB' sys/sys-mm.inp)" 1
"$C/cp2k-input.sh" -d sys -j qmmm -q LIG -Q 0 -n 2 2>/dev/null
check "QM/MM input"         "$([[ -s sys/sys-qmmm.inp ]] && echo yes)" yes
check "GEEP embedding"      "$(grep -c 'USE_GEEP_LIB' sys/sys-qmmm.inp)" 1
check "cubic QM cell"       "$(awk '/&QMMM/{q=1} q && /ABC/{print ($2==$3 && $3==$4) ? "yes" : "no"; exit}' sys/sys-qmmm.inp)" yes
check "basis for QM kinds"  "$(grep -c 'BASIS_SET DZVP-MOLOPT-SR-GTH' sys/sys-qmmm.inp | awk '{print ($1>0)?"yes":"no"}')" yes
CP2K=${CP2K_EXE:-cp2k.psmp}; command -v "$CP2K" >/dev/null || { echo "  cp2k not installed, run skipped"; finish; exit; }
(cd sys && "$CP2K" -i sys-mm.inp -o sys-mm.out > /dev/null 2>&1)
check "MM run finished"     "$(grep -c 'ENERGY| Total FORCE_EVAL ( FIST )' sys/sys-mm.out | awk '{print ($1>0)?"yes":"no"}')" yes
check "no abort"            "$(grep -ci abort sys/sys-mm.out)" 0
finish
