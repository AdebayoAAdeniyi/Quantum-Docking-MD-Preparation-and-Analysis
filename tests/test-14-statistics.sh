#!/usr/bin/env bash
# R statistics: PCA/correlation/regression on the real xTB screening table.
source "$(dirname "$0")/testlib.sh"
need Rscript
S=$QDM_ROOT/14-statistics
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT; cd "$W"
Rscript -e 'if (!requireNamespace("FactoMineR", quietly=TRUE)) quit(status=1)' 2>/dev/null ||
    { echo "  FactoMineR not installed, skipped"; exit 0; }
cp "$QDM_ROOT/examples/xtb-benzophenones/xtb-screen.csv" .
Rscript "$S/r-analysis/property-stats.R" -i xtb-screen.csv -o out > log.txt 2>&1
check "5 molecules read"        "$(awk '/molecules/{print $1}' log.txt)" 5
check "correlation matrix"      "$([[ -s out/correlation-matrix.csv ]] && echo yes)" yes
check "diagonal is 1"           "$(awk -F, 'NR==2{print $2}' out/correlation-matrix.csv)" 1
check "PCA variance sums to 100" "$(awk -F, 'NR>1{s+=$3} END{printf "%.1f", s}' out/pca-variance.csv)" 100.0 0.2
check "biplot written"          "$([[ -s out/pca-biplot.png ]] && echo yes)" yes
# add a response column and check the regression path
awk -F, 'NR==1{print $0",y"; next} {print $0","(NR*0.37)}' xtb-screen.csv > xy.csv
Rscript "$S/r-analysis/property-stats.R" -i xy.csv -y y -x "gap_eV,omega_eV,dipole_D" -o out2 > log2.txt 2>&1
check "regression ran"          "$([[ -s out2/regression.txt ]] && echo yes)" yes
check "correlations ranked"     "$(($(wc -l < out2/correlations-with-y.csv) - 1))" 3
check "R2 reported"             "$(grep -c 'R2 =' log2.txt)" 1
finish
