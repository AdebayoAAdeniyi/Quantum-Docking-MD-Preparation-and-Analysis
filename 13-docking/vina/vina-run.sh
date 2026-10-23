#!/usr/bin/env bash
# =============================================================================
#  vina-run.sh - dock every prepared ligand with AutoDock Vina
#  usage: vina-run.sh [-s vina-setup] [-e 16] [-n 9] [-c CPUS] [-S seed] [-q]
#    -e exhaustiveness (default 16)   -n poses per ligand (default 9)
#    -q write a job script (run-vina.sh) instead of docking now
#  Vina: Eberhardt, Santos-Martins, Tillack and Forli, J. Chem. Inf. Model. 2021,
#  61, 3891; Trott and Olson, J. Comput. Chem. 2010, 31, 455.
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"; source "$QDM_ROOT/lib/jobscript.sh"; load_cluster_conf
s=vina-setup; ex=16; np=9; cpu=${CORES:-4}; seed=42; queue=0
while getopts "s:e:n:c:S:qh" o; do case $o in s) s=$OPTARG ;; e) ex=$OPTARG ;; n) np=$OPTARG ;; c) cpu=$OPTARG ;; S) seed=$OPTARG ;; q) queue=1 ;; *) sed -n '2,9p' "$0"; exit 0 ;; esac; done
[[ -f "$s/receptor.pdbqt" && -f "$s/box.txt" ]] || die "run vina-prepare.sh first"
mkdir -p "$s/docked"
cmds=()
for l in "$s"/ligands/*.pdbqt; do
    n=$(basename "$l" .pdbqt)
    cmds+=("[[ -s docked/${n}_out.pdbqt ]] || vina --receptor receptor.pdbqt --ligand ligands/$n.pdbqt --config box.txt --exhaustiveness $ex --num_modes $np --cpu $cpu --seed $seed --out docked/${n}_out.pdbqt > docked/$n.log 2>&1")
done
((${#cmds[@]})) || die "no ligands in $s/ligands"
CORES=$cpu; JOB_THREADED=1
write_jobscript vina "$(cd "$s" && pwd)" VINA_SETUP "${cmds[@]}" > "$s/run-vina.sh"; chmod +x "$s/run-vina.sh"
if ((queue)); then info "$s/run-vina.sh written (${#cmds[@]} ligands)"
else (cd "$s" && SCHEDULER=local bash run-vina.sh > run-vina.log 2>&1) && info "docking finished: $s/docked" || die "vina failed ($s/run-vina.log)"; fi
