#!/usr/bin/env bash
# =============================================================================
#  aimall-run.sh - QTAIM analysis of many wavefunctions with AIMAll (no GUI)
#  usage: aimall-run.sh [-l list] [-p parent] [-s suffix] [-q]
#    looks for PARENT/NAME/NAME$SUFFIX.{wfx,fchk,wfn} (suffix default empty)
#    -q  write run-aimall.sh job scripts instead of running now
#  AIMAll (T. A. Keith, TK Gristmill Software) is commercial; set AIMQB_EXE in
#  config/cluster.conf or the environment (default: aimqb).
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"; source "$QDM_ROOT/lib/jobscript.sh"; load_cluster_conf
list=""; parent=gaussian-jobs; suf=""; queue=0
while getopts "l:p:s:qh" o; do case $o in l) list=$OPTARG ;; p) parent=$OPTARG ;; s) suf=$OPTARG ;; q) queue=1 ;; *) sed -n '2,9p' "$0"; exit 0 ;; esac; done
AIMQB=${AIMQB_EXE:-aimqb}; ((queue)) || require "$AIMQB"
if [[ -n "$list" ]]; then mapfile -t names < <(read_list "$list"); else mapfile -t names < <(find "$parent" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort); fi
for n in "${names[@]}"; do
    w=""; for e in wfx fchk wfn; do [[ -f "$parent/$n/$n$suf.$e" ]] && { w="$n$suf.$e"; break; }; done
    [[ -n "$w" ]] || { warn "$n: no wfx/fchk/wfn"; continue; }
    cmd="$AIMQB -nogui -nproc=\$NCORES $w > ${w%.*}.qblog 2>&1"
    write_jobscript "aim-$n" "$PWD/$parent/$n" AIMALL_SETUP "$cmd" > "$parent/$n/run-aimall.sh"; chmod +x "$parent/$n/run-aimall.sh"
    if ((queue)); then info "$n: run-aimall.sh written"; else (cd "$parent/$n" && SCHEDULER=local bash run-aimall.sh >/dev/null 2>&1) && info "$n: done" || warn "$n: AIMAll failed"; fi
done
