#!/usr/bin/env bash
# =============================================================================
#  rosetta-symm-dock.sh - build a cyclic oligomer under Cn symmetry
#
#  Writes the symmetry definition, then either scores the symmetric assembly or
#  docks it with SymDockProtocol. One subunit goes in; the symmetry definition
#  generates the rest, so only the rigid body degrees of freedom between
#  subunits are sampled.
#
#  usage: rosetta-symm-dock.sh -p SUBUNIT.pdb [options]
#    -p FILE    the single subunit to oligomerise
#    -N N       number of subunits (default: 3, a homotrimer)
#    -r R       initial subunit-to-axis distance in Angstrom (default: 40)
#    -s FILE    use this symmetry definition instead of generating one
#    -o DIR     output directory (default: symm_dock)
#    -n N       number of structures (default: 1)
#    -S         score only, no docking: check the symmetry definition first
#    -t         restore_talaris_behavior, for comparing against older results
#    -q         write the job script and stop, for the queue
#    -h         this help
#
#  Two things to know before using this.
#
#  The symmetry definition is geometry, and the numbers in it are not generic.
#  A definition written for one protein places that protein's subunits; reusing
#  it for a different one starts the search from a meaningless arrangement. With
#  -s you supply your own, which is what you want once you have tuned it. Rosetta
#  ships make_symmdef_file_denovo.py for generating these, and this script calls
#  it when it is available, falling back to a generated Cn definition otherwise.
#
#  Score only, with -S, is worth running first. It costs seconds and tells you
#  whether the definition builds the oligomer you expected. The original script
#  in this project only ever scored, despite being named for docking, so the
#  distinction is now an explicit flag rather than a surprise.
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../../lib/common.sh"
source "$QDM_ROOT/config/cluster.conf"
[[ -f "$QDM_ROOT/config/cluster.local.conf" ]] && source "$QDM_ROOT/config/cluster.local.conf"

pdb=""; nsub=3; radius=40; symdef=""; out="symm_dock"; nstruct=1
scoreonly=0; talaris=0; queue=0
while getopts "p:N:r:s:o:n:Stqh" o; do case $o in
    p) pdb=$OPTARG ;; N) nsub=$OPTARG ;; r) radius=$OPTARG ;; s) symdef=$OPTARG ;;
    o) out=$OPTARG ;; n) nstruct=$OPTARG ;; S) scoreonly=1 ;; t) talaris=1 ;;
    q) queue=1 ;;
    h) sed -n '2,/^# ====/p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) die "unknown option; -h for usage" ;;
esac; done

[[ -n "$pdb" ]] || die "a subunit PDB is required (-p); -h for usage"
[[ -f "$pdb" ]] || die "not found: $pdb"
[[ "$nsub" =~ ^[0-9]+$ ]] && (( nsub >= 2 )) || die "-N must be 2 or more"

# One chain only. A subunit that already contains several chains will be
# replicated whole, giving N times as many chains as intended.
chains=$(awk '/^(ATOM|HETATM)/{c[substr($0,22,1)]=1} END{print length(c)}' "$pdb")
(( chains <= 1 )) || die "$pdb has $chains chains; a symmetric subunit must be one chain
Split it first with prepare/rosetta-chains.sh"

mkdir -p "$out/input_files" "$out/output_files"
base=$(strip_ext "$pdb")
cp "$pdb" "$out/input_files/$base.pdb"

if [[ -n "$symdef" ]]; then
    [[ -f "$symdef" ]] || die "symmetry definition not found: $symdef"
    cp "$symdef" "$out/input_files/symm.symdef"
    info "using the supplied symmetry definition"
else
    gen=""
    for cand in "${ROSETTA_BIN:+$ROSETTA_BIN/../src/apps/public/symmetry/make_symmdef_file_denovo.py}" \
                "${ROSETTA_TOOLS:+$ROSETTA_TOOLS/../source/src/apps/public/symmetry/make_symmdef_file_denovo.py}" \
                "$(command -v make_symmdef_file_denovo.py 2>/dev/null || true)"; do
        [[ -n "$cand" && -f "$cand" ]] && { gen="$cand"; break; }
    done
    if [[ -n "$gen" ]]; then
        info "generating a C$nsub definition with Rosetta's make_symmdef_file_denovo.py"
        python "$gen" -symm_type cn -nsub "$nsub" -radius "$radius" \
            > "$out/input_files/symm.symdef" 2>/dev/null ||
        python "$gen" -symm_type cn -nsub "$nsub" \
            > "$out/input_files/symm.symdef" ||
        die "make_symmdef_file_denovo.py failed"
    else
        warn "make_symmdef_file_denovo.py not found; writing a plain C$nsub definition"
        warn "set ROSETTA_BIN or ROSETTA_TOOLS to use Rosetta's generator instead"
        {
            printf 'symmetry_name c%s\n' "$nsub"
            printf 'E = %s*VRT0_base' "$nsub"
            for ((i = 1; i < nsub; i++)); do printf ' + %s*(VRT0_base:VRT%s_base)' "$nsub" "$i"; done
            printf '\nanchor_residue COM\n'
            printf 'virtual_coordinates_start\n'
            awk -v n="$nsub" 'BEGIN{
                pi = atan2(0, -1)
                for (i = 0; i < n; i++) {
                    a = 2 * pi * i / n
                    cx = cos(a); cy = sin(a)
                    printf "xyz VRT%d  %.7f,%.7f,0.0000000  %.7f,%.7f,0.0000000  0.0000000,0.0000000,0.0000000\n", i, -cx, -cy, cy, -cx
                    printf "xyz VRT%d_base  %.7f,%.7f,0.0000000  %.7f,%.7f,0.0000000  0.0000000,0.0000000,0.0000000\n", i, -cx, -cy, cy, -cx
                }
                print "xyz VRT  -1.0000000,0.0000000,0.0000000  0.0000000,-1.0000000,0.0000000  0.0000000,0.0000000,0.0000000"
            }'
            printf 'virtual_coordinates_stop\n'
            for ((i = 0; i < nsub; i++)); do
                printf 'connect_virtual JUMP%s_to_com VRT%s VRT%s_base\n' "$i" "$i" "$i"
                printf 'connect_virtual JUMP%s_to_subunit VRT%s_base SUBUNIT\n' "$i" "$i"
            done
            for ((i = 0; i < nsub; i++)); do printf 'connect_virtual JUMP%s VRT VRT%s\n' "$i" "$i"; done
            printf 'set_dof JUMP0_to_com x(%s)\n' "-$radius"
            printf 'set_dof JUMP0_to_subunit angle_x(0) angle_y(0) angle_z(0)\n'
            printf 'set_jump_group JUMPGROUP2'
            for ((i = 0; i < nsub; i++)); do printf ' JUMP%s_to_com' "$i"; done; printf '\n'
            printf 'set_jump_group JUMPGROUP3'
            for ((i = 0; i < nsub; i++)); do printf ' JUMP%s_to_subunit' "$i"; done; printf '\n'
        } > "$out/input_files/symm.symdef"
    fi
fi

app=$( ((scoreonly)) && echo score || echo SymDock )
require_rosetta "$app"     # needed to run; the definition above is already written
exe=$(rosetta_exe "$app")

{
    cat <<EOF
# generated by rosetta-symm-dock.sh
-in:file:s input_files/$base.pdb
-ignore_unrecognized_res

# symmetry
-symmetry:symmetry_definition input_files/symm.symdef
-symmetry:initialize_rigid_body_dofs
-symmetry:symmetric_rmsd

# rotamers
-packing:ex1
-packing:ex2aro

-out:path:all output_files
-nstruct $nstruct
-out:file:fullatom
-overwrite
EOF
    ((talaris)) && printf -- '-restore_talaris_behavior\n'
    ((scoreonly)) || printf -- '-docking:dock_pert 3 8\n'
} > "$out/input_files/options"

cat > "$out/run.sh" <<EOF
#!/usr/bin/env bash
set -euo pipefail
${ROSETTA_SETUP:-}
cd "\$(dirname "\$0")"
$exe @input_files/options -out:suffix _symm

shopt -s nullglob
sc=(output_files/*.sc output_files/*.fasc)
if ((\${#sc[@]} == 0)); then
    echo "no score file was written; the run produced nothing" >&2
    exit 1
fi
"\$QDM_ROOT/15-rosetta-modelling/extract/rosetta-scores.sh" -o output_files/scores \\
    -k $( ((scoreonly)) && echo total_score || echo I_sc ) output_files
EOF
chmod +x "$out/run.sh"

info "wrote $out/input_files/symm.symdef and options"
info "application: $exe ($( ((scoreonly)) && echo 'scoring only' || echo 'symmetric docking' ))"
if ((queue)); then
    printf '\nSubmit it with:\n  bin/qdm-submit %s\n' "$out"
else
    QDM_ROOT="$QDM_ROOT" "$out/run.sh"
fi
