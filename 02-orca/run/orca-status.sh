#!/usr/bin/env bash
# =============================================================================
#  orca-status.sh - state of every ORCA step in a batch
#
#  usage: orca-status.sh [-l list] [-p parent]
#  For each NAME-STEP.inp: not-run / running / normal / error, whether the
#  geometry optimisation converged, SCF convergence and imaginary frequencies.
#  Restart hint: an unconverged optimisation leaves NAME-STEP.xyz with the last
#  geometry; copy it over the start geometry (or use * xyzfile with it) and
#  resubmit. Add "! MORead" with %moinp "NAME-STEP.gbw" to reuse the orbitals.
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"
AWKO="$QDM_ROOT/lib/awk/orca.awk"
list=""; parent=orca-jobs
while getopts "l:p:h" o; do case $o in l) list=$OPTARG ;; p) parent=$OPTARG ;; *) sed -n '2,11p' "$0"; exit 0 ;; esac; done
if [[ -n "$list" ]]; then mapfile -t names < <(read_list "$list")
else mapfile -t names < <(find "$parent" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort); fi
printf '%-28s %-28s %-9s %-6s %s\n' step status opt_conv n_imag note
for n in "${names[@]}"; do
    for inp in "$parent/$n"/*.inp; do
        [[ -f "$inp" ]] || continue
        o="${inp%.inp}.out"; s=$(basename "${inp%.inp}")
        if [[ ! -s "$o" ]]; then printf '%-28s %-28s\n' "$s" not-run; continue; fi
        IFS=, read -r _ st oc _ _ _ _ _ _ _ _ _ _ _ ni _ <<< "$(awk -v mode=summary -f "$AWKO" "$o")"
        note=""; [[ $st == incomplete ]] && { st=running; note=$(tail -n 1 "$o" | cut -c1-50); }
        want=0; grep -qi "OptTS\|NEB-TS" "$inp" && want=1
        [[ -n "$ni" && "$ni" != "$want" ]] && note="$ni imaginary frequencies (expected $want) $note"
        printf '%-28s %-28s %-9s %-6s %s\n' "$s" "$st" "$oc" "$ni" "$note"
    done
done
