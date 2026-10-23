#!/usr/bin/env bash
# =============================================================================
#  jobscript.sh - write a job script for local, PBS or SLURM execution
#
#  write_jobscript NAME WORKDIR SETUP_VAR COMMAND...   > run.sh
#     NAME       job name
#     WORKDIR    directory the job runs in
#     SETUP_VAR  name of the config variable holding environment lines
#                (GAUSSIAN_SETUP, GROMACS_SETUP, ...), or "" for none
#     set JOB_THREADED=1 before calling for shared memory programs (Gaussian)
#     so SLURM gets --cpus-per-task instead of MPI tasks
#     COMMAND    the command line(s), already quoted
#
#  The scheduler, cores, memory and time come from config/cluster.conf (or
#  config/cluster.local.conf if it exists, or the file in $QDM_CLUSTER_CONF).
# =============================================================================
load_cluster_conf() {
    local f
    for f in "${QDM_CLUSTER_CONF:-}" "$QDM_ROOT/config/cluster.local.conf" "$QDM_ROOT/config/cluster.conf"; do
        [[ -n "$f" && -f "$f" ]] && { source "$f"; return; }
    done
}

write_jobscript() {
    local name="$1" workdir="$2" setupvar="$3"; shift 3
    load_cluster_conf
    local setup=""; [[ -n "$setupvar" ]] && setup="${!setupvar}"
    echo "#!/bin/bash"
    case "${SCHEDULER:-local}" in
        pbs)
            echo "#PBS -N ${name:0:15}"
            [[ -n "$ACCOUNT" ]] && echo "#PBS -P $ACCOUNT"
            [[ -n "$QUEUE"   ]] && echo "#PBS -q $QUEUE"
            echo "#PBS -l select=${NODES:-1}:ncpus=${CORES:-1}:mpiprocs=${CORES:-1}${GPUS:+:ngpus=$GPUS}" | sed 's/:ngpus=0//'
            echo "#PBS -l walltime=${WALLTIME:-24:00:00}"
            echo "#PBS -o $name.pbs.out"
            echo "#PBS -e $name.pbs.err"
            [[ -n "$EMAIL" ]] && { echo "#PBS -m abe"; echo "#PBS -M $EMAIL"; }
            ;;
        slurm)
            echo "#SBATCH --job-name=$name"
            [[ -n "$ACCOUNT" ]] && echo "#SBATCH --account=$ACCOUNT"
            [[ -n "$QUEUE"   ]] && echo "#SBATCH --partition=$QUEUE"
            echo "#SBATCH --nodes=${NODES:-1}"
            if ((${JOB_THREADED:-0})); then        # shared memory programs (Gaussian, ORCA with OpenMP ...)
                echo "#SBATCH --ntasks=1"; echo "#SBATCH --cpus-per-task=${CORES:-1}"
            else
                echo "#SBATCH --ntasks-per-node=${CORES:-1}"
            fi
            echo "#SBATCH --mem=${MEM:-4G}"
            echo "#SBATCH --time=${WALLTIME:-24:00:00}"
            (( ${GPUS:-0} > 0 )) && echo "#SBATCH --gpus-per-node=$GPUS"
            echo "#SBATCH --output=$name.%j.out"
            [[ -n "$EMAIL" ]] && { echo "#SBATCH --mail-type=END,FAIL"; echo "#SBATCH --mail-user=$EMAIL"; }
            ;;
        local) echo "# runs on this machine: bash run.sh" ;;
        *) echo "unknown SCHEDULER=$SCHEDULER" >&2; return 1 ;;
    esac
    cat <<EOS
set -euo pipefail
export OMP_NUM_THREADS=\${OMP_NUM_THREADS:-1}
NCORES=${CORES:-1}
# run where the job was submitted, so job folders can be moved or copied
cd "\${PBS_O_WORKDIR:-\${SLURM_SUBMIT_DIR:-\$(cd "\$(dirname "\$0")" && pwd)}}"
${setup}
echo "job \$(basename "\$PWD") started \$(date)"
EOS
    printf '%s\n' "$@"
    echo 'echo "job finished $(date)"'
}

# submit_job DIR [SCRIPT]: submit DIR/run.sh with the configured scheduler
submit_job() {
    local dir="$1" script="${2:-run.sh}"
    load_cluster_conf
    ( cd "$dir" || exit 1
      case "${SCHEDULER:-local}" in
          pbs)   qsub   "$script" ;;
          slurm) sbatch "$script" ;;
          local) bash   "$script" > "${script%.sh}.log" 2>&1 && echo "finished $dir" || echo "FAILED $dir (see ${script%.sh}.log)" ;;
      esac )
}
