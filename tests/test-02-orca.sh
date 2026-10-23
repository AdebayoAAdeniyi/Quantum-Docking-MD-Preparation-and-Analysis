#!/usr/bin/env bash
# ORCA module: inputs for real molecules; parsers checked against public ORCA
# 4.2/5.0/6.0 outputs from cclib, plus two documented-layout fragments.
source "$(dirname "$0")/testlib.sh"
need obabel
O=$QDM_ROOT/02-orca; AWKO=$QDM_ROOT/lib/awk/orca.awk
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT; cd "$W"
for v in 5.0 6.0; do for f in dvb_ir.out dvb_td.out dvb_nmr.out dvb_raman.out dvb_gopt.out; do
    fetch "$CCLIB/ORCA/basicORCA$v/$f" "orca$v/$f" || exit 0; done; done
fetch "$CCLIB/ORCA/basicORCA6.0/dvb_sp_un_dft.out" orca6.0/dvb_sp_un_dft.out
fetch "$CCLIBDATA/ORCA/ORCA5.0/dvb_soc.log" orca5.0/dvb_soc.log
fetch "$CCLIBDATA/ORCA/ORCA4.2/casscf.out" orca4.2/casscf.out
fetch "$CCLIBDATA/ORCA/ORCA4.0/IrCl6_sp.out" orca4.0/IrCl6_sp.out

echo "== inputs"
cp "$QDM_ROOT"/examples/molecules/benzophenones/*.xyz .; printf '1-H\n2-OMe\n3-Oct\n4-Cl\n5-Br\n' > list
"$O/prepare/orca-input.sh" -l list -j opt,temps,td,soc,nmr -s smd:methanol -n 4 -M 16G 2>/dev/null
check "5 inputs per molecule"      "$(ls orca-jobs/3-Oct/*.inp | wc -l)" 5
check "maxcore = 75% of 16G / 4"   "$(grep -m1 maxcore orca-jobs/1-H/1-H-opt.inp | awk '{print $2}')" 3072
check "td starts from opt geometry" "$(tail -1 orca-jobs/1-H/1-H-td.inp)" "* xyzfile 0 1 1-H-opt.xyz"
check "SMD block"                  "$(grep -c 'smd true' orca-jobs/1-H/1-H-nmr.inp)" 1
"$O/prepare/orca-input.sh" -n 4 -o cas -j casscf -a 6,5 -r 4,3 -u 1,3 -S 1-H.xyz 2>/dev/null
check "NEVPT2 + SOC block"         "$(grep -c 'dosoc true' cas/1-H/1-H-cas3.inp)" 1
cp 4-Cl.xyz rx-r.xyz; cp 4-Cl.xyz rx-p.xyz
"$O/prepare/orca-input.sh" -n 4 -o neb -j neb rx 2>/dev/null
check "NEB-TS input"               "$(grep -c 'NEB_End_XYZFile "rx-p.xyz"' neb/rx/rx-neb.inp)" 1

echo "== parsers (ORCA 5 and 6 layouts)"
s() { awk -v mode="$1" -v name=t -f "$AWKO" "$DATA/$2"; }
check "G, ORCA 5"               "$(s summary orca5.0/dvb_ir.out | cut -d, -f10)" -381.91112705
check "G, ORCA 6"               "$(s summary orca6.0/dvb_ir.out | cut -d, -f10)" -381.91114546
check "opt converged"           "$(s summary orca6.0/dvb_gopt.out | cut -d, -f3)" yes
check "SCF failure flagged"     "$(s summary orca4.0/IrCl6_sp.out | cut -d, -f2)" error/SCF-not-converged
check "S1 f, ORCA 5"            "$(s uv orca5.0/dvb_td.out | head -1 | cut -d, -f6)" 0.005267213
check "S1 f, ORCA 6 (new table)" "$(s uv orca6.0/dvb_td.out | awk -F, '$3=="Singlet"' | head -1 | cut -d, -f6)" 0.005038390
check "S1 transitions"          "$(s uv orca6.0/dvb_td.out | awk -F, '$3=="Singlet"' | head -1 | cut -d, -f7)" '"33a->35a (56%); 34a->36a (43%)"'
check "IR intensity mode 8"     "$(s freq orca6.0/dvb_ir.out | awk -F, '$2==8{print $4}')" 0.36
check "Raman activity mode 7"   "$(s freq orca6.0/dvb_raman.out | awk -F, '$2==7{print $5}')" 2.788888
check "C1 shielding"            "$(s nmr orca6.0/dvb_nmr.out | head -1 | cut -d, -f4)" 114.398
check "Hirshfeld C1 (radical)"  "$(s charges orca6.0/dvb_sp_un_dft.out | head -1 | cut -d, -f6)" 0.069911
check "largest SOCME T9-S9"     "$(s socme orca5.0/dvb_soc.log | sort -t, -k4 -gr | head -1 | cut -d, -f4)" 12.18
check "SOC corrected state 1"   "$(s socuv orca5.0/dvb_soc.log | head -1 | cut -d, -f4)" 462.1
check "NEVPT2 S1 (eV)"          "$(s casscf orca4.2/casscf.out 2>/dev/null | awk -F, '$4==1{print $6}')" 3.505
check "NEB CI marked"           "$(awk -v mode=neb -f "$AWKO" "$TESTS_DIR/fixtures/orca-neb-fragment.out" | awk -F, '$6=="CI"{print $2}')" 6
check "fluorescence rate"       "$(awk -v mode=esd -f "$AWKO" "$TESTS_DIR/fixtures/orca-esd-fragment.out" | cut -d, -f3)" 1.688355e+06

echo "== batch extraction and ESD set-up"
mkdir -p jobs/dvb jobs/em
cp "$DATA/orca6.0/dvb_ir.out" jobs/dvb/dvb-opt.out; cp "$DATA/orca6.0/dvb_td.out" jobs/dvb/dvb-td.out
cp "$DATA/orca6.0/dvb_nmr.out" jobs/dvb/dvb-nmr.out
"$O/extract/orca-extract.sh" -k all -p jobs -o res -r "C=186.0,H=31.8" 2>/dev/null
check "summary rows"            "$(($(wc -l < res/summary.csv) - 1))" 3
check "shift column added"      "$(head -1 res/nmr.csv | awk -F, '{print $NF}')" shift_ppm
# mechanics only: a singlet and a doublet output stand in for the S0 and T1 jobs
cp "$DATA/orca6.0/dvb_ir.out" jobs/em/em-GS.out; cp "$DATA/orca6.0/dvb_sp_un_dft.out" jobs/em/em-TG.out
touch jobs/em/em-GS.hess jobs/em/em-TG.hess; awk -v mode=geometry -f "$AWKO" "$DATA/orca6.0/dvb_ir.out" > jobs/em/em-GS.xyz
echo em > emlist; "$O/prepare/orca-esd-input.sh" -p jobs -l emlist 2>/dev/null
dele=$(awk 'BEGIN{printf "%.1f", (-381.831980903245 - -382.055107107616) * 219474.6313632}')
check "DELE in cm-1"            "$(grep -m1 DELE jobs/em/em-phosp.inp | awk '{print $2}')" "$dele"
check "three T1 sublevels"      "$(grep -c 'IRoot' jobs/em/em-phosp.inp)" 3
need gnuplot
mkdir -p jobs/rx; cp "$TESTS_DIR/fixtures/orca-neb-fragment.out" jobs/rx/rx-neb.out
"$O/extract/orca-extract.sh" -k neb -p jobs -o res 2>/dev/null
"$O/plot/plot-neb.sh" -i res/neb.csv -o res/neb.png 2>/dev/null
check "NEB plot"                "$([[ -s res/neb.png ]] && echo yes)" yes
finish
