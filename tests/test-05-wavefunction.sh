#!/usr/bin/env bash
# Wavefunction analysis: QTAIM table logic on a layout fragment; Multiwfn and
# AIMAll are not installed in the test environment, so their drivers are only
# checked for syntax and option handling.
source "$(dirname "$0")/testlib.sh"
M=$QDM_ROOT/05-wavefunction-analysis
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT; cd "$W"
mkdir -p jobs/cplx; cp "$TESTS_DIR/fixtures/aimall-bcp-fragment.sum" jobs/cplx/cplx.sum
"$M/extract/qtaim-extract.sh" -p jobs -o res -e Fe 2>/dev/null
check "3 BCPs (ring CP ignored)"   "$(($(wc -l < res/bcp-all.csv) - 1))" 3
check "H = G + V"                  "$(awk -F, '$3=="C1"{print $10}' res/bcp-all.csv)" -0.8
check "covalent C-O"               "$(awk -F, '$3=="C1"{print $12}' res/bcp-all.csv)" shared
check "Fe-N intermediate"          "$(awk -F, '$3=="Fe3"{print $12}' res/bcp-metal.csv)" intermediate
check "weak contact filter"        "$(awk -F, 'NR>1{print $3"-"$4}' res/bcp-weak.csv)" H5-O6
for s in "$M"/run/*.sh "$M"/extract/*.sh; do bash -n "$s" && ok=yes || ok=no; check "syntax $(basename "$s")" "$ok" yes; done
finish
