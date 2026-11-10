#!/usr/bin/env bash
# =============================================================================
#  rosetta-cm-setup.sh - stage 1 of RosettaCM: alignments and membrane topology
#
#  Turns a multiple sequence alignment of the target against its templates into
#  one Grishin alignment file per template, which is the format RosettaCM's
#  partial_thread reads, and converts a predicted transmembrane topology into
#  the .span file the membrane score terms need.
#
#  usage: rosetta-cm-setup.sh -t TARGET -a ALIGNMENT [options]
#    -t TARGET     target name, matching the fasta basename (e.g. e2)
#    -a ALIGNMENT  alignment of target and templates, Clustal or two-column
#    -n NAME       name of the target *inside* the alignment (default: TARGET)
#    -f FASTA      target sequence (default: TARGET.fasta)
#    -O OCTOPUS    OCTOPUS or SPOCTOPUS topology file, for a membrane protein
#    -o DIR        output directory (default: 1_setup)
#    -h            this help
#
#  The alignment comes from a multiple aligner such as Clustal Omega, with the
#  target and every template in one file. Sequences must be upper case: Rosetta
#  treats lower case as a deletion marker, so mixed case silently truncates the
#  model.
#
#  Fragment files, TARGET_3.frags and TARGET_9.frags, are not produced here.
#  They come from the Robetta fragment server and are copied in if present.
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../../lib/common.sh"
[[ -f "$QDM_ROOT/config/cluster.local.conf" ]] && source "$QDM_ROOT/config/cluster.local.conf"

target=""; align=""; tname=""; fasta=""; octopus=""; out="1_setup"
while getopts "t:a:n:f:O:o:h" o; do case $o in
    t) target=$OPTARG ;; a) align=$OPTARG ;; n) tname=$OPTARG ;;
    f) fasta=$OPTARG ;; O) octopus=$OPTARG ;; o) out=$OPTARG ;;
    h) sed -n '2,/^# ====/p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) die "unknown option; -h for usage" ;;
esac; done

[[ -n "$target" ]] || die "a target name is required (-t); -h for usage"
[[ -n "$align"  ]] || die "an alignment file is required (-a); -h for usage"
[[ -f "$align"  ]] || die "alignment not found: $align"
tname="${tname:-$target}"
fasta="${fasta:-$target.fasta}"
[[ -f "$fasta" ]] || die "target fasta not found: $fasta"

mkdir -p "$out"
cp "$fasta" "$out/$target.fasta"

# Rosetta reads lower case as "delete this residue", so a sequence that came
# back from an aligner in mixed case would be silently truncated.
if grep -v '^>' "$fasta" | grep -q '[a-z]'; then
    warn "$fasta contains lower case; Rosetta reads that as a deletion"
    warn "uppercasing the copy in $out/"
    awk '/^>/{print; next}{print toupper($0)}' "$fasta" > "$out/$target.fasta"
fi

# Fragment libraries from the Robetta fragment server, if they are to hand.
for n in 3 9; do
    if [[ -f "${target}_${n}.frags" ]]; then cp "${target}_${n}.frags" "$out/"
    else warn "${target}_${n}.frags not found; get it from the Robetta fragment server"
    fi
done

# ---------------------------------------------------------------------------
# Membrane topology -> span
# ---------------------------------------------------------------------------
if [[ -n "$octopus" ]]; then
    [[ -f "$octopus" ]] || die "octopus file not found: $octopus"
    cp "$octopus" "$out/$target.octopus"
    # octopus2span.pl is a Rosetta tool by Bjorn Wallner and Vladimir
    # Yarov-Yarovoy. It is not redistributed here; point ROSETTA_TOOLS at a
    # Rosetta install. See docs/third-party-material.md.
    o2s=""
    for cand in "${ROSETTA_TOOLS:+$ROSETTA_TOOLS/membrane_tools/octopus2span.pl}" \
                "${ROSETTA_TOOLS:+$ROSETTA_TOOLS/octopus2span.pl}" \
                "$(command -v octopus2span.pl 2>/dev/null || true)"; do
        [[ -n "$cand" && -f "$cand" ]] && { o2s="$cand"; break; }
    done
    if [[ -n "$o2s" ]]; then
        perl "$o2s" "$out/$target.octopus" > "$out/$target.span"
        info "wrote $out/$target.span ($(grep -c . "$out/$target.span") lines)"
    else
        warn "octopus2span.pl not found; set ROSETTA_TOOLS in config/cluster.local.conf"
        warn "no .span file written, so the membrane score terms cannot be used"
    fi
fi

# ---------------------------------------------------------------------------
# Alignment -> one Grishin file per template
#
# The original script pulled each sequence out with "grep -i NAME", which
# matches substrings: a template called 2kz would also match 2kzq, and the
# wrong sequence would be threaded without any error. Names are compared as
# whole fields here instead.
# ---------------------------------------------------------------------------
# Accept Clustal output or a plain "NAME SEQUENCE" table. Clustal writes the
# alignment in blocks, so sequences are accumulated per name and joined.
awk '
    /^CLUSTAL|^MUSCLE|^MAFFT/ { next }                 # aligner banner
    /^[[:space:]]*$/          { next }                 # block separator
    /^[[:space:]]/            { next }                 # conservation line
    NF >= 2 {
        name = $1; seq = $2
        if (!(name in S)) order[++n] = name
        S[name] = S[name] seq
    }
    END { for (i = 1; i <= n; i++) printf "%s\t%s\n", order[i], S[order[i]] }
' "$align" > "$out/aligned.tsv"

grep -q . "$out/aligned.tsv" || die "no sequences parsed from $align; expected Clustal or 'NAME SEQUENCE' lines"

target_seq=$(awk -F'\t' -v n="$tname" '$1 == n {print $2}' "$out/aligned.tsv")
[[ -n "$target_seq" ]] || die "target '$tname' not found in $align
names present: $(awk -F'\t' '{printf "%s ", $1}' "$out/aligned.tsv")
use -n to give the name as it appears in the alignment"

# Every aligned sequence must be the same length, or the alignment is broken
# and the threading will be nonsense rather than merely poor.
if awk -F'\t' '{print length($2)}' "$out/aligned.tsv" | sort -u | wc -l | grep -qv '^1$'; then
    die "aligned sequences differ in length; $align is not a valid alignment"
fi

n_templates=0
while IFS=$'\t' read -r name seq; do
    [[ "$name" == "$tname" ]] && continue
    # Grishin format: a header naming the pair, a score line, then the two
    # gapped sequences each preceded by their start offset (0 based).
    {
        printf '## %s %s.pdb\n' "$tname" "$name"
        printf '#\n'
        printf 'scores from program: 0\n'
        printf '0 %s\n' "$target_seq"
        printf '0 %s\n' "$seq"
    } > "$out/${tname}_${name}.grishin"
    n_templates=$((n_templates + 1))
    info "wrote $out/${tname}_${name}.grishin"
done < "$out/aligned.tsv"

(( n_templates > 0 )) || die "no templates found in the alignment besides the target"
info "stage 1 complete: $n_templates template alignment(s) in $out/"
printf '\nNext: thread each template onto the target\n'
printf '  %s -t %s -o 2_threading -s %s/ TEMPLATE:CHAIN ...\n' \
       "15-rosetta-modelling/prepare/rosetta-cm-thread.sh" "$target" "$out"
