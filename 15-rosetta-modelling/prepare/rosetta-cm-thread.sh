#!/usr/bin/env bash
# =============================================================================
#  rosetta-cm-thread.sh - stage 2 of RosettaCM: thread templates onto the target
#
#  Runs partial_thread once per template, producing the partial models that the
#  hybridize stage recombines.
#
#  usage: rosetta-cm-thread.sh -t TARGET [options] TEMPLATE[:CHAIN] ...
#    -t TARGET     target name (e.g. e2)
#    -n NAME       target name inside the alignment (default: TARGET)
#    -s DIR        stage 1 output directory (default: 1_setup)
#    -d DIR        where the template PDB files are (default: .)
#    -o DIR        output directory (default: 2_threading)
#    -c            clean each template with Rosetta's clean_pdb.py first
#    -h            this help
#
#  Templates are named on the command line, optionally with the chain to use:
#      rosetta-cm-thread.sh -t e2 6urh:C 2kzq:A 4mwf:C
#  The original script carried this list inside the file, which meant editing
#  the script for every target.
#
#  On cleaning: clean_pdb.py is a Baker Lab and Meiler Lab script that ships
#  with Rosetta and is not redistributed here. With -c, ROSETTA_TOOLS must
#  point at a Rosetta install. Without it, the template PDB files are used as
#  given, which is what you want if you have already trimmed them by hand in a
#  viewer and renumbered them from 1. Rosetta needs each template numbered
#  continuously from 1 with a single chain; it will not tell you if that is
#  untrue, it will thread badly.
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../../lib/common.sh"
[[ -f "$QDM_ROOT/config/cluster.local.conf" ]] && source "$QDM_ROOT/config/cluster.local.conf"
[[ -n "${ROSETTA_SETUP:-}" ]] && eval "$ROSETTA_SETUP"

target=""; tname=""; setup="1_setup"; pdbdir="."; out="2_threading"; clean=0
while getopts "t:n:s:d:o:ch" o; do case $o in
    t) target=$OPTARG ;; n) tname=$OPTARG ;; s) setup=$OPTARG ;;
    d) pdbdir=$OPTARG ;; o) out=$OPTARG ;; c) clean=1 ;;
    h) sed -n '2,/^# ====/p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) die "unknown option; -h for usage" ;;
esac; done
shift $((OPTIND - 1))

[[ -n "$target" ]] || die "a target name is required (-t); -h for usage"
(($# > 0)) || die "name at least one template, e.g. 6urh:C; -h for usage"
tname="${tname:-$target}"
[[ -d "$setup" ]] || die "stage 1 directory not found: $setup (run rosetta-cm-setup.sh first)"
require_rosetta partial_thread
thread=$(rosetta_exe partial_thread)

mkdir -p "$out"
cp "$setup/$target.fasta" "$out/"
cp "$setup"/*.grishin "$out/" 2>/dev/null || die "no .grishin files in $setup"

threaded=0
for spec in "$@"; do
    code="${spec%%:*}"
    chain="${spec#*:}"; [[ "$chain" == "$spec" ]] && chain=""

    src=""
    for cand in "$pdbdir/$code.pdb" "$pdbdir/${code^^}.pdb" "$code.pdb"; do
        [[ -f "$cand" ]] && { src="$cand"; break; }
    done
    [[ -n "$src" ]] || { warn "template PDB not found for '$code' in $pdbdir; skipped"; continue; }

    grishin="$out/${tname}_${code}.grishin"
    [[ -f "$grishin" ]] || { warn "no alignment $grishin for template $code; skipped"; continue; }

    tpl="$out/$code.pdb"
    if ((clean)); then
        cp="${ROSETTA_TOOLS:+$ROSETTA_TOOLS/protein_tools/scripts/clean_pdb.py}"
        [[ -n "$cp" && -f "$cp" ]] || cp="$(command -v clean_pdb.py 2>/dev/null || true)"
        if [[ -n "$cp" && -f "$cp" ]]; then
            ( cd "$out" && cp "../$src" . 2>/dev/null || true
              python2 "$cp" "$(basename "$src")" "${chain:-A}" >/dev/null 2>&1 || \
              python  "$cp" "$(basename "$src")" "${chain:-A}" >/dev/null 2>&1 || true )
            # clean_pdb.py writes CODE_CHAIN.pdb with the code upper cased
            found=$(ls "$out/${code^^}_${chain:-A}.pdb" 2>/dev/null | head -1 || true)
            if [[ -n "$found" ]]; then mv "$found" "$tpl"
            else warn "clean_pdb.py produced nothing for $code; using the file as given"; cp "$src" "$tpl"; fi
        else
            warn "clean_pdb.py not found; set ROSETTA_TOOLS. Using $code as given"
            cp "$src" "$tpl"
        fi
    else
        cp "$src" "$tpl"
    fi

    info "threading $code onto $target"
    ( cd "$out" && "$thread" \
        -in:file:fasta "$target.fasta" \
        -in:file:alignment "${tname}_${code}.grishin" \
        -in:file:template_pdb "$code.pdb" ) || die "partial_thread failed on $code"

    # partial_thread appends .pdb to the template name
    if [[ -f "$out/$code.pdb.pdb" ]]; then
        mv "$out/$code.pdb.pdb" "$out/${tname}_on_$code.pdb"
        info "wrote $out/${tname}_on_$code.pdb"
        threaded=$((threaded + 1))
    else
        warn "partial_thread wrote no model for $code"
    fi
done

(( threaded > 0 )) || die "no templates were threaded; nothing for stage 3 to recombine"
info "stage 2 complete: $threaded threaded model(s) in $out/"
printf '\nNext: recombine them\n'
printf '  15-rosetta-modelling/run/rosetta-cm-hybridize.sh -t %s -s %s -n 10000\n' "$target" "$out"
