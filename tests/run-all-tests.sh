#!/usr/bin/env bash
# =============================================================================
#  run-all-tests.sh - run every module test, skipping what is not installed
#
#  usage: tests/run-all-tests.sh [-s] [NUMBER ...]
#    -s   also run the slow parts (QDM_SLOW=1): longer QM and MD runs
#    NUMBER  run only these modules, e.g. "01 04 13"
#  Tests that need a program which is not installed report "skipped" and do not
#  fail the run. Tests that need the internet download small public reference
#  outputs once into tests/data (cached afterwards).
# =============================================================================
cd "$(dirname "$0")"
while getopts "sh" o; do case $o in s) export QDM_SLOW=1 ;; *) sed -n '2,11p' "$0"; exit 0 ;; esac; done
shift $((OPTIND-1))
sel=("$@")
pass=0; fail=0; failed=()
for t in test-*.sh; do
    n=$(echo "$t" | sed 's/test-\([0-9]*\).*/\1/')
    if ((${#sel[@]})); then printf '%s\n' "${sel[@]}" | grep -qx "$n" || continue; fi
    printf '\n\033[1m=== %s\033[0m\n' "$t"
    if bash "$t"; then pass=$((pass + 1)); else fail=$((fail + 1)); failed+=("$t"); fi
done
printf '\n========================================\n'
printf 'modules passed: %d, failed: %d\n' "$pass" "$fail"
((fail)) && printf 'failed: %s\n' "${failed[*]}"
exit $((fail > 0))
