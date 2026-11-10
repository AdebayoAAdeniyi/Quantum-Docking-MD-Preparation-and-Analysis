#!/usr/bin/env bash
# =============================================================================
#  rosetta-chains.sh - relabel and renumber the chains of an assembled PDB
#
#  Takes a PDB whose chains are separated by TER records and gives each one a
#  distinct chain ID, optionally renumbering every chain from 1. Assemblies
#  built up in a viewer, or concatenated from separate models, routinely arrive
#  with every chain labelled A or with no chain ID at all, and Rosetta's
#  partner syntax (AB_C) cannot address them until that is fixed.
#
#  usage: rosetta-chains.sh -p FILE.pdb [options]
#    -p FILE    input PDB, chains separated by TER
#    -o FILE    output PDB (default: INPUT_chains.pdb)
#    -c IDS     chain IDs to assign in order (default: ABCDEFGHIJ)
#    -r         renumber each chain from 1
#    -h         this help
#
#  On -r: renumbering per chain is what you want for symmetric docking, where
#  each subunit should look identical, and for PDB-numbered restraints, where
#  163A and 163B should mean the equivalent residue in two copies. It is wrong
#  if you need to keep the numbering of a reference structure, so it is not the
#  default.
#
#  The original version of this script did the renumbering by passing each chain
#  through tleap, which rebuilds the residues and can change atom names and
#  protonation on the way. This does it by arithmetic on the residue sequence
#  number field, so nothing else about the file changes. It also wrote two
#  outputs, one renumbered and one not, with similar names; this writes one.
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../../lib/common.sh"

pdb=""; out=""; ids="ABCDEFGHIJ"; renumber=0
while getopts "p:o:c:rh" o; do case $o in
    p) pdb=$OPTARG ;; o) out=$OPTARG ;; c) ids=$OPTARG ;; r) renumber=1 ;;
    h) sed -n '2,/^# ====/p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) die "unknown option; -h for usage" ;;
esac; done

[[ -n "$pdb" ]] || die "an input PDB is required (-p); -h for usage"
[[ -f "$pdb" ]] || die "not found: $pdb"
out="${out:-$(strip_ext "$pdb")_chains.pdb}"

n_ter=$(grep -c '^TER' "$pdb" || true)
n_atoms=$(grep -c '^\(ATOM\|HETATM\)' "$pdb" || true)
(( n_atoms > 0 )) || die "$pdb contains no ATOM or HETATM records"
(( n_ter > 0 )) || die "$pdb has no TER records, so the chain boundaries are unknown
Add a TER after each chain, including the last one, or split the file first."

info "$n_atoms atoms, $n_ter TER record(s), so $n_ter chain(s) expected"
(( ${#ids} >= n_ter )) || die "only ${#ids} chain ID(s) given for $n_ter chain(s); use -c"

awk -v ids="$ids" -v renum="$renumber" '
    BEGIN { ci = 1 }
    /^(ATOM|HETATM)/ {
        id = substr(ids, ci, 1)
        if (renum) {
            resseq = substr($0, 23, 4)
            if (resseq != last_seq) { last_seq = resseq; newnum[ci]++ }
            printf "%s%s%4d%s\n", substr($0, 1, 21), id, newnum[ci], substr($0, 27)
        } else {
            printf "%s%s%s\n", substr($0, 1, 21), id, substr($0, 23)
        }
        next
    }
    /^TER/ { print "TER"; ci++; last_seq = ""; next }
    /^END/ { next }
    { print }
    END { print "END" }
' "$pdb" > "$out"

got=$(awk '/^(ATOM|HETATM)/{c[substr($0,22,1)]=1} END{for(k in c) printf "%s", k}' "$out" | fold -w1 | sort | tr -d '\n')
kept=$(grep -c '^\(ATOM\|HETATM\)' "$out")
(( kept == n_atoms )) || die "atom count changed, $n_atoms in and $kept out; not writing over anything"

info "wrote $out"
info "chains: $got$( ((renumber)) && echo ', each renumbered from 1' )"
awk '/^(ATOM|HETATM)/{c=substr($0,22,1); n[c]++;
      r=substr($0,23,4)+0; if(!(c in lo) || r<lo[c]) lo[c]=r; if(r>hi[c]) hi[c]=r}
     END{for(k in n) printf "  chain %s: %5d atoms, residues %d to %d\n", k, n[k], lo[k], hi[k]}' \
    "$out" | sort
