#!/usr/bin/env bash
# GAMESS-US / Firefly: mixed basis + ECP inputs from the Basis Set Exchange, NEDA
# fragment detection, parsers checked on public cclib outputs and the NBO tutorial.
source "$(dirname "$0")/testlib.sh"
need obabel; need bse
M=$QDM_ROOT/03-gamess-firefly; AWKM=$QDM_ROOT/lib/awk/gamess.awk
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT; cd "$W"
for f in basicGAMESS-US2017/dvb_ir.out basicGAMESS-US2017/dvb_gopt_a.out basicGAMESS-US2017/dvb_un_sp.out \
         basicGAMESS-US2017/MoOCl4-sp.out basicGAMESS-US2017/MoOCl4-sp.inp basicGAMESS-US2017/nh3_ts_ir.out basicFirefly8.0/dvb_ir.out; do
    fetch "$CCLIB/GAMESS/$f" "gamess/$f" || exit 0
done
echo "== inputs"
# MoOCl4 geometry from the cclib reference; our Mo basis/ECP must equal the one in that working input
awk -v mode=geometry -f "$AWKM" "$DATA/gamess/basicGAMESS-US2017/MoOCl4-sp.out" > MoOCl4.xyz
echo "MoOCl4 -2 1" > chargemult.txt
"$M/prepare/gamess-input.sh" -j energy -n 2 -M 4G MoOCl4.xyz 2>/dev/null
i=gamess-jobs/MoOCl4/MoOCl4-energy.inp
check "PP=READ for the ECP"           "$(grep -c 'PP=READ' $i)" 1
check "Mo first exponent = reference" "$(awk '/^ Mo /{getline; getline; print $2; exit}' $i)" "$(awk '/^ Mo /{getline; getline; print $2; exit}' "$DATA/gamess/basicGAMESS-US2017/MoOCl4-sp.inp")"
check "Mo ECP first term = reference" "$(awk '/MO-ECP GEN/{getline; getline; print $1; exit}' $i)" -0.0469492
check "light atoms get ECP NONE"      "$(grep -c 'CL-ECP NONE' $i)" 4
obabel -:"O.O" -oxyz -O wd.xyz --gen3d 2>/dev/null
"$M/prepare/gamess-input.sh" -P firefly -j neda -f "1,3,4/2,5,6" -n 2 -M 4G -o ff wd.xyz 2>/dev/null
check "NEDA fragments (Firefly)"      "$(grep NEDA ff/wd/wd-neda.inp | xargs)" "NAO NEDA (1,3,4) (2,5,6) END"
check "Firefly P2P group"             "$(grep -c P2P ff/wd/wd-neda.inp)" 1
echo "== parsers"
s() { awk -v mode="$1" -f "$AWKM" "$DATA/gamess/$2"; }
check "E final (GAMESS-US)"  "$(s summary basicGAMESS-US2017/dvb_ir.out | cut -d, -f8)" -382.0506353313
check "ZPE (GAMESS-US)"      "$(s summary basicGAMESS-US2017/dvb_ir.out | cut -d, -f9)" 0.177076
check "HOMO (eV)"            "$(s summary basicGAMESS-US2017/dvb_ir.out | cut -d, -f12)" -4.0654
check "Firefly CRLF output"  "$(s summary basicFirefly8.0/dvb_ir.out | cut -d, -f2,3)" "Firefly,normal"
check "opt converged"        "$(s summary basicGAMESS-US2017/dvb_gopt_a.out | cut -d, -f4)" converged
check "SCF failure flagged"  "$(s summary basicGAMESS-US2017/dvb_un_sp.out | cut -d, -f3)" normal/SCF-not-converged
check "TS: one imaginary"    "$(s summary basicGAMESS-US2017/nh3_ts_ir.out | cut -d, -f16)" 1
check "TS mode"              "$(s freq basicGAMESS-US2017/nh3_ts_ir.out | head -1 | cut -d, -f3)" -825.18
check "dianion charge"       "$(s summary basicGAMESS-US2017/MoOCl4-sp.out | cut -d, -f5)" -2
check "NEDA total (NBO tutorial)" "$(awk -v mode=neda -f "$AWKM" "$TESTS_DIR/fixtures/nbo-neda-fragment.out" | awk -F, '$2=="Total_interaction"{print $3}')" -6.22
check "NEDA EL"              "$(awk -v mode=neda -f "$AWKM" "$TESTS_DIR/fixtures/nbo-neda-fragment.out" | awk -F, '$2 ~ /^EL/{print $3}')" -14.25
echo "== batch"
mkdir -p jobs/wd; cp "$TESTS_DIR/fixtures/nbo-neda-fragment.out" jobs/wd/wd-neda.log
"$M/extract/gamess-extract.sh" -k neda -p jobs -o res 2>/dev/null
check "wide NEDA table"      "$(awk -F, 'NR==2{print NF}' res/neda-wide.csv)" 14
finish
