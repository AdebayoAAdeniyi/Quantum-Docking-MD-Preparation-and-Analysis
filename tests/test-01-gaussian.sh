#!/usr/bin/env bash
# Gaussian module: input generation for real molecules, and extraction checked
# against public Gaussian 16 / 09 outputs from the cclib project
# (https://github.com/cclib/cclib, BSD-3 licence; downloaded, not redistributed).
source "$(dirname "$0")/testlib.sh"
need obabel
G=$QDM_ROOT/01-gaussian; AWKG=$QDM_ROOT/lib/awk/gaussian.awk
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT; cd "$W"

echo "== input generation (benzophenones + a Mo complex with an ECP)"
cp "$QDM_ROOT"/examples/molecules/benzophenones/*.xyz .
for f in dvb_ir.out dvb_td.out dvb_nmr.log Trp_polar.log dvb_un_sp.log Mo4OCl4-sp.log dvb_gopt.out; do
    fetch "$CCLIB/Gaussian/basicGaussian16/$f" "g16/$f" || exit 0
done
for f in dvb_hirshfeld.out water_neutral_nbo_opt.out water_cation_nbo_opt.out water_neutral_nbo_opt_freq.out; do
    fetch "$CCLIBDATA/Gaussian/Gaussian16/$f" "g16/$f" || exit 0
done
fetch "$CCLIBDATA/Gaussian/Gaussian09/Dahlgren_TS.log" g09/Dahlgren_TS.log
fetch "$CCLIBDATA/Gaussian/Gaussian09/dvb_gopt_unconverged.log" g09/dvb_gopt_unconverged.log
fetch "$CCLIBDATA/Gaussian/Gaussian03/tms_nmr.log" g03/tms_nmr.log
cp "$DATA/g16/Mo4OCl4-sp.log" MoOCl4.log
printf '1-H\n2-OMe\n3-Oct\n4-Cl\n5-Br\nMoOCl4\n' > list; echo "MoOCl4 -1 2" > chargemult.txt
"$G/prepare/gaussian-input.sh" -l list -j opt,pop,chelpg,td,nmr,polar -s smd:acetonitrile -F -n 4 -M 8G 2>/dev/null
check "6 job folders"                "$(ls -d gaussian-jobs/*/ | wc -l)" 6
check "5 Link1 steps in 1-H.com"     "$(grep -c -- '--Link1--' gaussian-jobs/1-H/1-H.com)" 5
check "Gen basis only in step 1 (Mo)" "$(grep -c '^Mo 0' gaussian-jobs/MoOCl4/MoOCl4.com)" 2
check "ChkBasis in later steps (Mo)" "$(grep -c 'ChkBasis' gaussian-jobs/MoOCl4/MoOCl4.com)" 5
check "Mo listed first"              "$(awk 'NR==9{print $1}' gaussian-jobs/MoOCl4/MoOCl4.com)" Mo
check "anion charge in fukui input"  "$(grep -A2 'N+1 (anion)' gaussian-jobs/1-H/1-H-fukui.com | tail -1)" "-1 2"
check "impossible multiplicity stops" "$(echo '1-H 0 2' > bad.txt; "$G/prepare/gaussian-input.sh" -c bad.txt -o x 1-H.xyz >/dev/null 2>&1; echo $?)" 1
"$G/prepare/gaussian-followup.sh" -k temps -T "298.15 310.15" -l list 2>/dev/null
check "Freq=ReadFC follow-up"        "$(grep -c 'Freq=ReadFC' gaussian-jobs/5-Br/5-Br-temps.com)" 2

echo "== extraction from real Gaussian outputs"
s() { awk -v mode="$1" -v name=t -f "$AWKG" "$DATA/$2"; }
check "G (dvb_ir)"           "$(s summary g16/dvb_ir.out | cut -d, -f12)" -382.164915
check "ZPE (dvb_ir)"         "$(s summary g16/dvb_ir.out | cut -d, -f9)"  0.177132
check "lowest freq (dvb_ir)" "$(s summary g16/dvb_ir.out | cut -d, -f20)" 53.1981
check "S1 energy (dvb_td)"   "$(s uv g16/dvb_td.out | head -1 | cut -d, -f4)" 5.3351
check "S1 transitions"       "$(s uv g16/dvb_td.out | head -1 | cut -d, -f7)" '"34->36 (46%); 35->36 (19%); 35->37 (34%)"'
check "C1 shielding"         "$(s nmr g16/dvb_nmr.log | head -1 | cut -d, -f4)" 114.0880
check "J(C3-H7) total"       "$(awk -v mode=jcoupling -v pairs=3-7 -f "$AWKG" "$DATA/g16/dvb_nmr.log" | cut -d, -f10)" 133.415
check "alpha iso (Trp)"      "$(s polar g16/Trp_polar.log | cut -d, -f2)" 74.1242
check "beta parallel (Trp)"  "$(s polar g16/Trp_polar.log | cut -d, -f9)" 51.5144
check "Fermi contact C1"     "$(s hyperfine g16/dvb_un_sp.log | head -1 | cut -d, -f5)" 28.58292
check "Hirshfeld C1"         "$(s charges g16/dvb_hirshfeld.out | head -1 | cut -d, -f7)" 0.001436
check "CM5 C1"               "$(s charges g16/dvb_hirshfeld.out | head -1 | cut -d, -f8)" -0.009905
check "NPA O, cation (total density block)" "$(s charges g16/water_cation_nbo_opt.out | head -1 | cut -d, -f6)" 0.12011
check "CCSD(T) energy"       "$(s summary g09/Dahlgren_TS.log | cut -d, -f7)" -434.37573219
check "unconverged opt flagged" "$(s summary g09/dvb_gopt_unconverged.log | cut -d, -f2)" error
check "Mo complex charge"    "$(s summary g16/Mo4OCl4-sp.log | cut -d, -f3)" -2

echo "== batch extraction, Link1 splitting, Fukui table and plots"
mkdir -p jobs/dvb jobs/water
cat "$DATA"/g16/{dvb_gopt.out,dvb_ir.out,dvb_td.out,dvb_nmr.log} > jobs/dvb/dvb.log
cp "$DATA/g16/water_neutral_nbo_opt_freq.out" jobs/water/water.log
# mechanics only: real N and N-1 outputs, N reused as N+1 (no anion file available)
cat "$DATA"/g16/{water_neutral_nbo_opt.out,water_neutral_nbo_opt.out,water_cation_nbo_opt.out} > jobs/water/water-fukui.log
printf 'dvb\nwater\n' > l2
"$G/extract/gaussian-extract.sh" -k all -l l2 -p jobs -o res -R "$DATA/g03/tms_nmr.log" 2>/dev/null
check "4 steps found in dvb.log" "$(grep -c '^dvb,' res/summary.csv)" 4
check "IP of water (eV)"         "$(awk -F, 'NR==2{print $2}' res/fukui-global.csv)" 8.0589
check "NMR shift uses TMS ref"   "$(awk -F, 'NR==2{print $8}' res/nmr.csv)" 189.4802
need gnuplot
"$QDM_ROOT/bin/qdm-plot-spectra" -t ir -i res/freq.csv -s 0.967 -o res/ir.png 2>/dev/null
"$QDM_ROOT/bin/qdm-plot-spectra" -t uv -i res/uv.csv -x 150:300 -o res/uv.png 2>/dev/null
check "IR plot written" "$([[ -s res/ir.png ]] && echo yes)" yes
check "UV plot written" "$([[ -s res/uv.png ]] && echo yes)" yes
finish
