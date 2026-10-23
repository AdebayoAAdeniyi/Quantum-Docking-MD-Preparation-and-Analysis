#!/usr/bin/env bash
# NWChem and xTB: real calculations run here (both programs are free).
#   xTB: GFN2 optimisation + Hessian, vertical IP/EA and Fukui indices for two benzophenones
#   NWChem: B3LYP/6-31G* optimisation + frequencies of water at two temperatures
# Set QDM_SLOW=1 to add the HCN -> HNC transition state search and the water dimer
# counterpoise job (about two minutes on two cores).
source "$(dirname "$0")/testlib.sh"
need obabel
X=$QDM_ROOT/04-nwchem-xtb
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT; cd "$W"
if command -v xtb >/dev/null; then
    echo "== xTB"
    cp "$QDM_ROOT"/examples/molecules/benzophenones/{1-H,4-Cl}.xyz .
    "$X/run/xtb-run.sh" -s water -n 2 1-H.xyz 4-Cl.xyz 2>/dev/null
    "$X/extract/xtb-extract.sh" 2>/dev/null
    check "two molecules screened"    "$(($(wc -l < xtb-results/xtb-screen.csv) - 1))" 2
    check "no imaginary modes"        "$(awk -F, 'NR>1{s+=$8}END{print s}' xtb-results/xtb-screen.csv)" 0
    check "IP > EA"                   "$(awk -F, 'NR>1 && $9>$10{n++}END{print n}' xtb-results/xtb-screen.csv)" 2
    check "omega = mu^2/(2 eta)"      "$(awk -F, 'NR==2{printf "%.3f", $11*$11/(2*$12) - $13}' xtb-results/xtb-screen.csv)" 0.000 0.001
    check "Fukui rows = atoms of 1-H" "$(grep -c '^1-H,' xtb-results/xtb-fukui.csv)" "$(head -1 1-H.xyz)"
else echo "  xtb not installed, part skipped"; fi
if command -v nwchem >/dev/null; then
    echo "== NWChem"
    obabel -:"O" -oxyz -O water.xyz --gen3d 2>/dev/null
    "$X/prepare/nwchem-input.sh" -j opt -b "6-31g*" -T "298.15 310" -n 1 -M 2G water.xyz 2>/dev/null
    (cd nwchem-jobs/water && mkdir -p perm && nwchem water-opt.nw > water-opt.log 2>&1)
    "$X/extract/nwchem-extract.sh" -k all 2>/dev/null
    check "water optimised"        "$(awk -F, 'NR==2{print $3,$4}' nwchem-results/summary.csv)" "normal yes"
    check "E(B3LYP/6-31G*) water"  "$(awk -F, 'NR==2{print $7}' nwchem-results/summary.csv)" -76.40702 0.00002
    check "3 vibrations, none imaginary" "$(awk -F, 'NR>1 && $4>0{n++}END{print n}' nwchem-results/freq.csv)" 3
    check "two temperatures"       "$(($(wc -l < nwchem-results/thermo.csv) - 1))" 2
    check "G = H - TS (298 K)"     "$(awk -F, 'NR==2{printf "%.5f", $9 - ($6 - $3*$7/1000/627.509474)}' nwchem-results/thermo.csv)" 0.00000 0.00001
    if [[ -n "${QDM_SLOW:-}" ]]; then
        printf '3\nTS guess\nC 0 0 0\nN 0 0 1.19\nH 1.1 0 0.65\n' > hcn-ts.xyz
        "$X/prepare/nwchem-input.sh" -j ts -b "6-31g*" -T 298.15 -n 1 -M 2G hcn 2>/dev/null
        (cd nwchem-jobs/hcn && mkdir -p perm && nwchem hcn-ts.nw > hcn-ts.log 2>&1)
        check "HCN/HNC TS: one imaginary mode" "$(awk -v mode=summary -f "$QDM_ROOT/lib/awk/nwchem.awk" nwchem-jobs/hcn/hcn-ts.log | cut -d, -f13)" 1
        "$X/run/nwchem-displace.sh" nwchem-jobs/hcn/hcn-ts.log disp 2>/dev/null
        check "displaced geometries" "$(ls disp/hcn-plus.xyz disp/hcn-minus.xyz | wc -l)" 2
    fi
else echo "  nwchem not installed, part skipped"; fi
finish
