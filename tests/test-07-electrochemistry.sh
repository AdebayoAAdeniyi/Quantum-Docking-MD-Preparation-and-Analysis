#!/usr/bin/env bash
# Electrochemistry: CV and UV processing on simulated data with analytically known
# answers, and the Born-Haber arithmetic against a hand calculation.
source "$(dirname "$0")/testlib.sh"
E=$QDM_ROOT/07-electrochemistry
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT; cd "$W"
echo "== cyclic voltammetry (simulated)"
awk -v type=fc -v e12=450 -f "$TESTS_DIR/fixtures/simulate-cv.awk" > FC.TXT
awk -v type=red -v e0=-1700 -v s=30 -f "$TESTS_DIR/fixtures/simulate-cv.awk" > 1-H.TXT
"$E/experimental/cv-process.sh" -F FC.TXT -d red -g 3.6 -O cv 1-H.TXT 2>/dev/null
check "Fc E1/2 = 0.450 V"     "$(awk -F, 'NR==2{print $4}' cv/ferrocene.csv)" 0.450 0.003
check "Fc dEp = 70 mV"        "$(awk -F, 'NR==2{print $5}' cv/ferrocene.csv)" 70 3
# the wave is i = -ilim/(1+exp((E-e0)/s)); its tangent-intersection onset is e0 + 2s = -1640 mV,
# which is -1640 - 450 = -2090 mV = -2.090 V vs Fc+/Fc
check "onset vs Fc = -2.090 V" "$(awk -F, 'NR==2{print $3}' cv/onsets.csv)" -2.090 0.015
check "E_LUMO = -(onset+4.8)"  "$(awk -F, 'NR==2{print $6}' cv/onsets.csv)" -2.71 0.02
check "E_HOMO = LUMO - Eg"     "$(awk -F, 'NR==2{print $7}' cv/onsets.csv)" -6.31 0.02
echo "== UV-Vis (simulated two-band spectrum)"
awk 'BEGIN{print "nm Abs"; for (l=600;l>=250;l-=0.5) printf "%.1f %.6f\n", l, exp(-((l-350)/20)^2) + 0.6*exp(-((l-280)/15)^2)}' > uv.txt
"$E/experimental/uv-onset.sh" -O uv uv.txt 2>/dev/null
check "lambda_max = 350 nm"   "$(awk -F, 'NR==2{print $2}' uv/uv-onsets.csv)" 350.0 0.6
check "Eg = 1239.84/onset"    "$(awk -F, 'NR==2{printf "%.3f", 1239.84/$4 - $5}' uv/uv-onsets.csv)" 0.000 0.002
echo "== Born-Haber cycle"
# invented but exactly reproducible energies: G_gas, E_gas, E_solv in Hartree
cat > en.csv <<'EOT'
name,G_gas,E_gas,E_solv
OX,-100.000000,-100.020000,-100.060000
RED,-100.100000,-100.120000,-100.220000
EOT
echo "test OX RED 1" > rx.txt
"$E/computational/redox-pka.sh" -k redox -t en.csv -r rx.txt -R SHE -o . 2>/dev/null
# dG_solv(OX) = (-100.06 + 100.02)*627.509474 + 1.894 = -25.100 - ... ; computed here the same way
want=$(awk 'BEGIN{ K=627.509474; R=1.987204259e-3; T=298.15; std=R*T*log(0.0820573661*T)
    gsO=(-100.060000 - -100.020000)*K + std; gsR=(-100.220000 - -100.120000)*K + std
    GO=-100.000000*K + gsO; GR=-100.100000*K + gsR
    printf "%.3f", -(GR-GO)/23.06054783 - 4.28 }')
check "E vs SHE"              "$(awk -F, 'NR==2{print $10}' redox.csv)" "$want" 0.001
check "1 atm -> 1 M included" "$(awk -F, 'NR==2{printf "%.2f", $6 - ((-100.060000 - -100.020000)*627.509474)}' redox.csv)" 1.89 0.01
finish
