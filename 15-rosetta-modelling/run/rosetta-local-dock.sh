#!/usr/bin/env bash
# =============================================================================
#  rosetta-local-dock.sh - local refinement and local docking of an assembly
#
#  For a complex you have already placed, by symmetry, by hand in a viewer, or
#  from a homologous structure. Each interface is first refined to remove
#  clashes, then perturbed and re-docked locally, and the best pose is carried
#  into the next interface.
#
#  usage: rosetta-local-dock.sh -p COMPLEX.pdb -i PARTNERS [-i PARTNERS ...] [options]
#    -p FILE      the assembled complex
#    -i PARTNERS  interface in Rosetta partner syntax, e.g. AB_C. Repeat for
#                 each interface to treat, in order.
#    -o DIR       output directory (default: local_dock)
#    -n N         poses per docking stage (default: 10)
#    -t "R T"     perturbation, degrees and Angstrom (default: "3 8")
#    -R           refine only, no perturbation and no re-docking
#    -F           skip the final full-complex relax
#    -q           write the job script and stop, for the queue
#    -h           this help
#
#  Example, a heterotrimer treated one interface at a time:
#    rosetta-local-dock.sh -p e1_trimer.pdb -i AB_C -i AC_B -n 10
#
#  On -n: the original script used -nstruct 1 for the refinement and 2 for the
#  docking, then picked the best of two. Local docking with a 3 degree and 8
#  Angstrom perturbation is a stochastic search, so two poses is not a search,
#  it is two samples. Ten is a reasonable floor and hundreds is normal.
#
#  On ordering: each interface starts from the best pose of the previous one, so
#  the order you list them in matters and the result is not symmetric. For a
#  genuinely symmetric oligomer use rosetta-symm-dock.sh instead, which
#  constrains all the interfaces to stay identical.
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../../lib/common.sh"
source "$QDM_ROOT/config/cluster.conf"
[[ -f "$QDM_ROOT/config/cluster.local.conf" ]] && source "$QDM_ROOT/config/cluster.local.conf"

pdb=""; out="local_dock"; nstruct=10; pert="3 8"; refineonly=0; norelax=0; queue=0
partners=()
while getopts "p:i:o:n:t:RFqh" o; do case $o in
    p) pdb=$OPTARG ;; i) partners+=("$OPTARG") ;; o) out=$OPTARG ;;
    n) nstruct=$OPTARG ;; t) pert=$OPTARG ;; R) refineonly=1 ;; F) norelax=1 ;;
    q) queue=1 ;;
    h) sed -n '2,/^# ====/p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) die "unknown option; -h for usage" ;;
esac; done

[[ -n "$pdb" ]] || die "a complex PDB is required (-p); -h for usage"
[[ -f "$pdb" ]] || die "not found: $pdb"
(( ${#partners[@]} )) || die "name at least one interface with -i, e.g. -i AB_C"

# Check every chain named in the partner strings is actually in the file. A
# typo here is not an error in Rosetta, it just docks nothing.
present=$(awk '/^(ATOM|HETATM)/{c[substr($0,22,1)]=1} END{for(k in c) printf "%s", k}' "$pdb")
for spec in "${partners[@]}"; do
    [[ "$spec" == *_* ]] || die "interface '$spec' is not in partner syntax, e.g. AB_C"
    for ((i = 0; i < ${#spec}; i++)); do
        ch="${spec:i:1}"
        [[ "$ch" == "_" ]] && continue
        [[ "$present" == *"$ch"* ]] || die "interface '$spec' names chain $ch, which is not in $pdb
chains present: $present"
    done
done
info "chains in $pdb: $present; interfaces: ${partners[*]}"

mkdir -p "$out/input_files" "$out/output_files"
base=$(strip_ext "$pdb")
cp "$pdb" "$out/input_files/$base.pdb"

# Needed to run, not to check the chains and write the script.
require_rosetta docking_protocol
((norelax)) || require_rosetta relax
docker=$(rosetta_exe docking_protocol)
{
    printf '#!/usr/bin/env bash\nset -euo pipefail\n%s\ncd "$(dirname "$0")"\n' "${ROSETTA_SETUP:-}"
    printf 'current=input_files/%s.pdb\n\n' "$base"
    printf 'pick_best() {   # pick_best SCOREFILE TERM -> best decoy name\n'
    printf '    "$QDM_ROOT/15-rosetta-modelling/extract/rosetta-scores.sh" \\\n'
    printf '        -o output_files/scores_$3 -k "$2" -n 1 "$1" >/dev/null\n'
    printf '    tail -n +2 output_files/scores_$3/scores-sorted.csv | head -1 | cut -d, -f1\n'
    printf '}\n\n'

    stage=0
    for spec in "${partners[@]}"; do
        stage=$((stage + 1))
        tag="refine${stage}"
        cat <<EOF
# ---- interface $spec: local refinement, to clear clashes ----
cat > flags_$tag <<FLAGS
-in:file:s \$current
-nstruct 1
-partners $spec
-docking_local_refine
-use_input_sc
-ex1
-ex2aro
-overwrite
-out:file:fullatom
-out:path:all output_files
-out:suffix _$tag
FLAGS
$docker @flags_$tag
best=\$(pick_best output_files/score_$tag.sc total_score $tag)
[[ -n "\$best" ]] || { echo "refinement of $spec produced no scored pose" >&2; exit 1; }
current=output_files/\$best.pdb
echo "after $tag: \$current"

EOF
    done

    if ((! refineonly)); then
        stage=0
        for spec in "${partners[@]}"; do
            stage=$((stage + 1))
            tag="dock${stage}"
            cat <<EOF
# ---- interface $spec: local docking, perturbation $pert ----
cat > flags_$tag <<FLAGS
-in:file:s \$current
-nstruct $nstruct
-partners $spec
-dock_pert $pert
-ex1
-ex2aro
-overwrite
-out:file:fullatom
-out:path:all output_files
-out:suffix _$tag
FLAGS
$docker @flags_$tag
# Ranked on I_sc, the interface energy, not on total_score. A pose can have the
# better total score while making a worse interface, which is the opposite of
# what a docking run is being asked.
best=\$(pick_best output_files/score_$tag.sc I_sc $tag)
[[ -n "\$best" ]] || { echo "docking of $spec produced no scored pose" >&2; exit 1; }
current=output_files/\$best.pdb
echo "after $tag: \$current"

EOF
        done
    fi

    if ((! norelax)); then
        relaxer=$(rosetta_exe relax)
        cat <<EOF
# ---- final relax of the whole complex ----
cat > flags_final_relax <<FLAGS
-in:file:s \$current
-nstruct 1
-relax:dualspace
-relax:minimize_bond_angles
-set_weights cart_bonded 0.5 pro_close 0
-default_max_cycles 200
-overwrite
-out:file:fullatom
-out:pdb
-out:path:all output_files
-out:suffix _final_relax
FLAGS
$relaxer @flags_final_relax
EOF
    fi

    printf '\ncp "$current" output_files/%s_selected.pdb\n' "$base"
    printf 'echo\necho "selected pose: output_files/%s_selected.pdb"\n' "$base"
} > "$out/run.sh"
chmod +x "$out/run.sh"

info "wrote $out/run.sh: ${#partners[@]} interface(s), $nstruct pose(s) per docking stage"
if ((queue)); then
    printf '\nSubmit it with:\n  bin/qdm-submit %s\n' "$out"
else
    QDM_ROOT="$QDM_ROOT" "$out/run.sh"
fi
