#!/usr/bin/env bash
# =============================================================================
#  gaussian-status.sh - check a batch of Gaussian jobs and prepare restarts
#
#  For each molecule prints: finished steps / total steps, final status, the
#  number of imaginary frequencies and the Gaussian error (if any).
#
#  usage: gaussian-status.sh [-l list] [-p parent] [-r]
#    -l LIST     molecules (default: every folder in parent)
#    -p DIR      job folder (default gaussian-jobs)
#    -r          for failed or unfinished jobs, save the last geometry as
#                restart/NAME.xyz so they can be resubmitted with
#                gaussian-input.sh -d restart -l restart/list
#  Typical messages and what to do are listed in docs/troubleshooting.md.
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"
AWKG="$QDM_ROOT/lib/awk/gaussian.awk"
list=""; parent=gaussian-jobs; restart=0
while getopts "l:p:rh" o; do
    case $o in l) list=$OPTARG ;; p) parent=$OPTARG ;; r) restart=1 ;; *) sed -n '2,15p' "$0"; exit 0 ;; esac
done
if [[ -n "$list" ]]; then mapfile -t names < <(read_list "$list")
else mapfile -t names < <(find "$parent" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort); fi
((restart)) && { mkdir -p restart; : > restart/list; }
printf '%-22s %-7s %-11s %-6s %s\n' molecule steps status n_imag message
for n in "${names[@]}"; do
    com="$parent/$n/$n.com"; log="$parent/$n/$n.log"
    total=0; [[ -f "$com" ]] && total=$(( $(grep -ci -- '--link1--' "$com"; true) + 1 ))
    if [[ ! -f "$log" ]]; then printf '%-22s %-7s %-11s %-6s %s\n' "$n" "0/$total" "not-run" "" ""; continue; fi
    done_=$(grep -c "Normal termination" "$log"; true)
    # Gaussian adds a frequency step itself after "Opt Freq", count it too
    auto=$(grep -ciE "^#.*(opt.*freq|freq.*opt)" "$com" 2>/dev/null; true)
    total=$((total + ${auto:-0})); ((total)) || total=1
    if grep -q "Error termination" "$log"; then
        status=error
        msg=$(grep -m1 -oE "Optimization stopped|Number of steps exceeded|Convergence failure|SCF has not converged|Erroneous write|galloc: could not allocate memory|Problem with the distance matrix|Wrong number of Negative eigenvalues|FormBX had a problem|Error in internal coordinate system|Linear angle in Bend|The combination of multiplicity [0-9]+ and +[0-9]+ electrons is impossible|End of file in ZSymb|Atoms too close|Small interatomic distances" "$log")
        [[ -z "$msg" ]] && msg="stopped in $(grep -m1 -oE "l[0-9]+\.exe" <(grep "Error termination" "$log"))"
    elif ((done_ >= total)); then status=done; msg=""
    else status=running; msg=$(tail -n 1 "$log" | sed 's/^ *//' | cut -c1-70); fi
    nim=$(awk -v mode=summary -f "$AWKG" "$log" | cut -d, -f19)
    want=0; grep -qiE "opt[^ ]*[(,=](ts|qst)" "$com" 2>/dev/null && want=1      # a TS needs exactly one
    [[ -n "$nim" && "$nim" != "$want" ]] && msg="$nim imaginary frequencies, expected $want. $msg"
    printf '%-22s %-7s %-11s %-6s %s\n' "$n" "$done_/$total" "$status" "$nim" "$msg"
    if ((restart)) && [[ $status != done || ( -n "$nim" && "$nim" != "$want" ) ]]; then
        awk -v mode=geometry -v name="$n" -f "$AWKG" "$log" > "restart/$n.xyz"
        echo "$n" >> restart/list
    fi
done
((restart)) && [[ -s restart/list ]] && info "restart geometries in restart/ ($(wc -l < restart/list) molecules)"
true
