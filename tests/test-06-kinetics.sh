#!/usr/bin/env bash
# Kinetics from the real NWChem B3LYP/6-31G* results stored in examples/kinetics-hcn-hnc
source "$(dirname "$0")/testlib.sh"
E=$QDM_ROOT/examples/kinetics-hcn-hnc
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT; cd "$W"
"$QDM_ROOT/06-kinetics/kinetics.sh" -r "$E/reactions.txt" -s "$E/summary.csv" -t "$E/thermo.csv" -w -o k 2>/dev/null
b() { awk -F, -v c="$1" 'NR==1{for(i=1;i<=NF;i++) if($i==c) k=i} NR==2{print $k}' k/barriers.csv; }
check "dE# (kcal/mol)"     "$(b 'dE#_kcal')" 50.13
check "dG# (kcal/mol)"     "$(b 'dG#_kcal')" 45.44
check "dE_r (kcal/mol)"    "$(b 'dE_r_kcal')" 16.59
check "kJ = 4.184 kcal"    "$(b 'dG#_kJ')" 190.13
check "Wigner kappa"       "$(awk -F, 'NR==2{print $4}' k/rates.csv)" 2.2775
# Eyring by hand: kappa kB T/h exp(-dG#/RT)
want=$(awk 'BEGIN{printf "%.4e", 2.2775*1.380649e-23*298.15/6.62607015e-34*exp(-45.44*4184/(8.314462618*298.15))}')
check "Eyring rate"        "$(awk -F, 'NR==2{print $5}' k/rates.csv)" "$want" 2e-23
check "4 temperatures"     "$(($(wc -l < k/rates-vs-T.csv) - 1))" 4
check "Arrhenius R2"       "$(awk -F, 'NR==2{print $5}' k/arrhenius.csv)" 1.00000 0.0001
need gnuplot
"$QDM_ROOT/06-kinetics/plot-profile.sh" -i k/profile.dat -o k/profile.png 2>/dev/null
check "profile plot" "$([[ -s k/profile.png ]] && echo yes)" yes
finish
