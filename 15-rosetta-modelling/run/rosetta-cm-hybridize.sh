#!/usr/bin/env bash
# =============================================================================
#  rosetta-cm-hybridize.sh - stage 3 of RosettaCM: recombine the threaded models
#
#  Writes the RosettaScripts XML and the options file, then runs the Hybridize
#  mover, which samples recombinations of the threaded templates and closes the
#  gaps between them with fragment insertion.
#
#  usage: rosetta-cm-hybridize.sh -t TARGET [options]
#    -t TARGET     target name (e.g. e2)
#    -s DIR        stage 2 output directory (default: 2_threading)
#    -S DIR        stage 1 output directory, for fasta/frags/span (default: 1_setup)
#    -o DIR        output directory (default: 3_hybridize)
#    -n N          number of models to generate (default: 1000)
#    -M            membrane protein: use the membrane weight sets and the span file
#    -w DIR        weights directory (default: the module's wts/)
#    -c FILE       extra constraint file applied in all three stages
#    -q            write the job script and stop, for the queue
#    -h            this help
#
#  nstruct is a real decision, not a formality. RosettaCM explores by
#  recombination, so a few hundred models rarely finds the best arrangement of
#  several templates. The published runs behind this module used 10000 to 15000.
#  The default here is 1000 so that a first run finishes and can be inspected.
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../../lib/common.sh"
source "$QDM_ROOT/config/cluster.conf"
[[ -f "$QDM_ROOT/config/cluster.local.conf" ]] && source "$QDM_ROOT/config/cluster.local.conf"

target=""; thr="2_threading"; setup="1_setup"; out="3_hybridize"
nstruct=1000; membrane=0; wts=""; cst=""; queue=0
while getopts "t:s:S:o:n:Mw:c:qh" o; do case $o in
    t) target=$OPTARG ;; s) thr=$OPTARG ;; S) setup=$OPTARG ;; o) out=$OPTARG ;;
    n) nstruct=$OPTARG ;; M) membrane=1 ;; w) wts=$OPTARG ;; c) cst=$OPTARG ;;
    q) queue=1 ;;
    h) sed -n '2,/^# ====/p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) die "unknown option; -h for usage" ;;
esac; done

[[ -n "$target" ]] || die "a target name is required (-t); -h for usage"
[[ -d "$thr"    ]] || die "stage 2 directory not found: $thr"
wts="${wts:-$QDM_ROOT/15-rosetta-modelling/wts}"
[[ -d "$wts" ]] || die "weights directory not found: $wts"
[[ "$nstruct" =~ ^[0-9]+$ ]] && (( nstruct > 0 )) || die "-n must be a positive integer"

mkdir -p "$out"
cp "$setup/$target.fasta" "$out/" 2>/dev/null || die "$setup/$target.fasta not found"
for n in 3 9; do
    cp "$setup/${target}_${n}.frags" "$out/" 2>/dev/null || \
        die "${target}_${n}.frags not found in $setup; Hybridize needs both fragment files"
done

models=("$thr"/*_on_*.pdb)
[[ -e "${models[0]}" ]] || die "no threaded models (*_on_*.pdb) in $thr"
cp "${models[@]}" "$out/"

if ((membrane)); then
    cp "$setup/$target.span" "$out/" 2>/dev/null || \
        die "-M needs $setup/$target.span; run rosetta-cm-setup.sh with -O"
    s1=stage1_membrane.wts; s2=stage2_membrane.wts; s3=stage3_rlx_membrane.wts
else
    s1=stage1.wts; s2=stage2.wts; s3=stage3.wts
fi
for w in "$s1" "$s2" "$s3"; do
    cp "$wts/$w" "$out/" 2>/dev/null || die "weights file not found: $wts/$w"
done
[[ -n "$cst" ]] && { cp "$cst" "$out/extra.cst" || die "constraint file not found: $cst"; }

# ---------------------------------------------------------------------------
# The XML. Templates are listed one per threaded model; cst_file="AUTO" tells
# Hybridize to derive per-template coordinate constraints from the template
# itself, which is what keeps the well-modelled regions in place.
# ---------------------------------------------------------------------------
{
    cat <<EOF
<ROSETTASCRIPTS>
    <SCOREFXNS>
        <ScoreFunction name="stage1" weights="$s1" symmetric="0">
            <Reweight scoretype="atom_pair_constraint" weight="1"/>
        </ScoreFunction>
        <ScoreFunction name="stage2" weights="$s2" symmetric="0">
            <Reweight scoretype="atom_pair_constraint" weight="0.5"/>
        </ScoreFunction>
        <ScoreFunction name="fullatom" weights="$s3" symmetric="0">
            <Reweight scoretype="atom_pair_constraint" weight="0.5"/>
        </ScoreFunction>
    </SCOREFXNS>
    <MOVERS>
        <Hybridize name="hybridize" stage1_scorefxn="stage1" stage2_scorefxn="stage2"
                   fa_scorefxn="fullatom" batch="1" stage1_increase_cycles="1.0"
                   stage2_increase_cycles="1.0" linmin_only="1">
            <Fragments three_mers="${target}_3.frags" nine_mers="${target}_9.frags"/>
EOF
    for m in "${models[@]}"; do
        printf '            <Template pdb="%s" cst_file="AUTO" weight="1.000" />\n' "$(basename "$m")"
    done
    cat <<'EOF'
        </Hybridize>
    </MOVERS>
    <PROTOCOLS>
        <Add mover="hybridize"/>
    </PROTOCOLS>
</ROSETTASCRIPTS>
EOF
} > "$out/rosetta_cm.xml"

# ---------------------------------------------------------------------------
# The options file
# ---------------------------------------------------------------------------
{
    cat <<EOF
# RosettaCM hybridize, generated by rosetta-cm-hybridize.sh
-in:file:fasta $target.fasta
-parser:protocol rosetta_cm.xml
-nstruct $nstruct

# relax settings used inside the full atom stage
-relax:minimize_bond_angles
-relax:minimize_bond_lengths
-relax:jump_move true
-relax:min_type lbfgs_armijo_nonmonotone
-default_max_cycles 200
-score:weights $s3
-use_bicubic_interpolation
-hybridize:stage1_probability 1.0

# trim unused residue type patches to keep the memory footprint down
-chemical:exclude_patches LowerDNA UpperDNA Cterm_amidation SpecialRotamer VirtualBB ShoveBB VirtualDNAPhosphate VirtualNTerm CTermConnect sc_orbitals pro_hydroxylated_case1 pro_hydroxylated_case2 ser_phosphorylated thr_phosphorylated tyr_phosphorylated tyr_sulfated lys_dimethylated lys_monomethylated lys_trimethylated lys_acetylated glu_carboxylated cys_acetylated tyr_diiodinated N_acetylated C_methylamidated MethylatedProteinCterm

-out:file:silent_struct_type binary
-overwrite
EOF
    ((membrane)) && cat <<EOF

# membrane scoring
-membrane
-in:file:spanfile $target.span
-membrane:no_interpolate_Mpair
-membrane:Menv_penalties
EOF
    [[ -n "$cst" ]] && printf '\n-constraints:cst_file extra.cst\n-constraints:cst_fa_file extra.cst\n'
} > "$out/rosetta_cm.options"

info "wrote $out/rosetta_cm.xml and $out/rosetta_cm.options"
info "templates: ${#models[@]}, nstruct: $nstruct, membrane: $((membrane))"

# ---------------------------------------------------------------------------
# The run script. An MPI build is used when one is configured, because nstruct
# in the thousands is the whole cost of this stage and a serial binary leaves
# the rest of the node idle. The original PBS script asked for 24 MPI ranks and
# then called the serial binary, so 23 cores did nothing.
# ---------------------------------------------------------------------------
exe=$(rosetta_exe rosetta_scripts) || die "rosetta_scripts not found; set ROSETTA_BIN and ROSETTA_SUFFIX"
mpi_exe=$(rosetta_exe rosetta_scripts mpi 2>/dev/null || true)
if [[ -n "$mpi_exe" && "${CORES:-1}" -gt 1 ]]; then
    runline="mpirun -np ${CORES} $mpi_exe @rosetta_cm.options"
    info "using the MPI build with ${CORES} ranks"
else
    runline="$exe @rosetta_cm.options"
    [[ "${CORES:-1}" -gt 1 ]] && warn "no MPI build found; running serially on one core.
Set ROSETTA_MPI_SUFFIX in config/cluster.local.conf to use all ${CORES}."
fi

cat > "$out/run.sh" <<EOF
#!/usr/bin/env bash
set -euo pipefail
${ROSETTA_SETUP:-}
cd "\$(dirname "\$0")"
$runline
mkdir -p output_structures
# The original script ran "mv S_*.pdb output_structures/." unconditionally, so
# a crashed run ended with a confusing "cannot stat S_*.pdb" and no indication
# that the real problem was upstream.
shopt -s nullglob
models=(S_*.pdb)
if ((\${#models[@]})); then
    mv "\${models[@]}" output_structures/
    echo "moved \${#models[@]} model(s) into output_structures/"
else
    echo "no S_*.pdb models were produced" >&2
    echo "check the log above: this stage usually dies on memory, on a missing" >&2
    echo "fragment file, or on a template that is not numbered from 1." >&2
    exit 1
fi
EOF
chmod +x "$out/run.sh"

if ((queue)); then
    info "job script ready: $out/run.sh"
    printf '\nSubmit it with:\n  bin/qdm-submit %s\n' "$out"
else
    info "running hybridize now (use -q to write the script and submit instead)"
    "$out/run.sh"
    printf '\nNext: relax and rank the models\n'
    printf '  15-rosetta-modelling/run/rosetta-relax.sh -t %s -s %s/output_structures%s\n' \
           "$target" "$out" "$( ((membrane)) && echo ' -M' )"
fi
