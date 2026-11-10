#!/usr/bin/env bash
# Rosetta modelling: everything up to launching Rosetta itself.
#
# Rosetta is free for non-commercial users but is not installed here, so what
# is tested is alignment preparation, score parsing, restraint checking and
# chain handling. Those are where the defects described in
# 15-rosetta-modelling/README.md lived: the ranking that read the wrong column,
# the restraints that were silently dropped, and the substring template match.
source "$(dirname "$0")/testlib.sh"
R=$QDM_ROOT/15-rosetta-modelling
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT; cd "$W"

# ---------------------------------------------------------------------------
echo "== rosetta-cm-setup.sh: Grishin alignments"
# A Clustal-style alignment in blocks, with two template names where one is a
# prefix of the other. "grep -i 2kz" matches both, which is the original bug.
cat > t.aln <<'EOF'
CLUSTAL O(1.2.4) multiple sequence alignment

E2_target   MKVLATSAL
2kz         MKVLA-SAL
2kzq        MKV-ATSAL

E2_target   WGCFRE
2kz         WGCF-E
2kzq        WG-FRE
EOF
printf '>E2_target\nMKVLATSALWGCFRE\n' > e2.fasta
"$R/prepare/rosetta-cm-setup.sh" -t e2 -a t.aln -n E2_target -o setup >/dev/null 2>&1
check "one grishin per template"   "$(ls setup/*.grishin 2>/dev/null | wc -l)" 2
check "target not aligned to itself" "$(ls setup/E2_target_E2_target.grishin 2>/dev/null | wc -l)" 0
# The decisive check: 2kz must carry its own sequence, not 2kzq's.
check "substring name not confused" \
    "$(tail -1 setup/E2_target_2kz.grishin | awk '{print $2}')" "MKVLA-SALWGCF-E"
check "longer name keeps its own"  \
    "$(tail -1 setup/E2_target_2kzq.grishin | awk '{print $2}')" "MKV-ATSALWG-FRE"
check "target sequence on line 4"  \
    "$(sed -n '4p' setup/E2_target_2kz.grishin | awk '{print $2}')" "MKVLATSALWGCFRE"
check "grishin header names the pair" \
    "$(head -1 setup/E2_target_2kz.grishin)" "## E2_target 2kz.pdb"
check "blocks were joined, not truncated" \
    "$(tail -1 setup/E2_target_2kz.grishin | awk '{print length($2)}')" 15

echo "== rosetta-cm-setup.sh: rejects bad input"
cat > ragged.aln <<'EOF'
E2_target   MKVLATSAL
2kz         MKVLA
EOF
"$R/prepare/rosetta-cm-setup.sh" -t e2 -a ragged.aln -n E2_target -o bad1 >/dev/null 2>&1
check "unequal sequence lengths rejected" "$?" 1
"$R/prepare/rosetta-cm-setup.sh" -t e2 -a t.aln -n NOT_THERE -o bad2 >/dev/null 2>&1
check "missing target name rejected"      "$?" 1
# Lower case is a deletion marker to Rosetta, so it must not pass through.
printf '>E2_target\nmkvlatsalwgcfre\n' > lower.fasta
"$R/prepare/rosetta-cm-setup.sh" -t lower -a t.aln -n E2_target -f lower.fasta -o lc >/dev/null 2>&1
check "lower case sequence uppercased"    \
    "$(grep -v '^>' lc/lower.fasta | tr -d '\n')" "MKVLATSALWGCFRE"

# ---------------------------------------------------------------------------
echo "== rosetta-scores.sh: classic SCORE: layout"
cat > a.sc <<'EOF'
SEQUENCE: AAAA
SCORE: total_score rms fa_rep I_sc description
SCORE: -120.5 1.2 55.1 -5.2 model_0001
SCORE: -135.8 0.9 52.3 -8.7 model_0002
SCORE: -99.1 2.4 60.0 2.1 model_0003
EOF
"$R/extract/rosetta-scores.sh" -o s1 -k total_score a.sc >/dev/null 2>&1
check "header row not treated as a model" "$(( $(wc -l < s1/scores.csv) - 1 ))" 3
check "best by total_score"      "$(sed -n '2p' s1/scores-sorted.csv | cut -d, -f1)" model_0002
"$R/extract/rosetta-scores.sh" -o s2 -k I_sc a.sc >/dev/null 2>&1
check "best by I_sc"             "$(sed -n '2p' s2/scores-sorted.csv | cut -d, -f1)" model_0002
check "worst by I_sc is last"    "$(tail -1 s2/scores-sorted.csv | cut -d, -f1)" model_0003
"$R/extract/rosetta-scores.sh" -o s3 -k nonexistent_term a.sc >/dev/null 2>&1
check "unknown term rejected"    "$?" 1

echo "== rosetta-scores.sh: column position must not matter"
# Same data, one score term removed, so I_sc is no longer column 5. The old
# "sort -nk5" picks model_0001 here; ranking by name must still pick 0002.
cat > c.sc <<'EOF'
SCORE: total_score fa_rep I_sc description
SCORE: -120.5 55.1 -5.2 model_0001
SCORE: -135.8 52.3 -8.7 model_0002
EOF
check "old positional sort is wrong here" \
    "$(sort -nk5 c.sc | head -1 | awk '{print $NF}')" model_0001
"$R/extract/rosetta-scores.sh" -o s4 -k I_sc c.sc >/dev/null 2>&1
check "lookup by name is right"  "$(sed -n '2p' s4/scores-sorted.csv | cut -d, -f1)" model_0002

echo "== rosetta-scores.sh: JSON .fasc layout"
cat > b.fasc <<'EOF'
{"decoy": "relax_0001", "total_score": -140.22, "fa_rep": 50.1}
{"decoy": "relax_0002", "total_score": -131.70, "fa_rep": 51.9}
EOF
"$R/extract/rosetta-scores.sh" -o s5 -k total_score b.fasc >/dev/null 2>&1
check "JSON fasc parsed"         "$(( $(wc -l < s5/scores.csv) - 1 ))" 2
check "JSON fasc ranked"         "$(sed -n '2p' s5/scores-sorted.csv | cut -d, -f1)" relax_0001
check "descending sort reverses" \
    "$("$R/extract/rosetta-scores.sh" -o s6 -k total_score -r b.fasc >/dev/null 2>&1; sed -n '2p' s6/scores-sorted.csv | cut -d, -f1)" relax_0002

echo "== rosetta-scores.sh: pooled files with different columns"
"$R/extract/rosetta-scores.sh" -o s7 -k total_score a.sc b.fasc >/dev/null 2>&1
check "pooling is reported"      "$([[ -s s7/incomplete-terms.txt ]] && echo yes)" yes

# ---------------------------------------------------------------------------
echo "== rosetta-chains.sh"
python3 -c "
lines=[]
for ch in range(3):
    for r in range(1,4):
        for an in (' N  ',' CA ',' C  '):
            lines.append('ATOM  %5d %4s ALA A%4d    %8.3f%8.3f%8.3f  1.00  0.00           C' % (len(lines)+1, an, 100+r+ch*10, r*1.5, ch*2.0, 0.5))
    lines.append('TER')
open('tri.pdb','w').write('\n'.join(lines)+'\nEND\n')
"
"$R/prepare/rosetta-chains.sh" -p tri.pdb -r -o fixed.pdb >/dev/null 2>&1
check "atom count preserved"     "$(grep -c '^ATOM' fixed.pdb)" 27
check "three chain IDs assigned" "$(awk '/^ATOM/{c[substr($0,22,1)]=1} END{print length(c)}' fixed.pdb)" 3
check "each chain renumbered from 1" \
    "$(awk '/^ATOM/{print substr($0,22,1), substr($0,23,4)+0}' fixed.pdb | sort -u | awk '$2==1' | wc -l)" 3
"$R/prepare/rosetta-chains.sh" -p tri.pdb -c AB -o few.pdb >/dev/null 2>&1
check "too few chain IDs rejected" "$?" 1
grep -v '^TER' tri.pdb > noter.pdb
"$R/prepare/rosetta-chains.sh" -p noter.pdb -o x.pdb >/dev/null 2>&1
check "missing TER records rejected" "$?" 1

# ---------------------------------------------------------------------------
echo "== rosetta-constrained-relax.sh: restraint checking"
cat > good.cst <<'EOF'
2A:CA   2B:CA   1.5  6.5
3A:N    3C:N    1.5  6.5
EOF
"$R/run/rosetta-constrained-relax.sh" -p fixed.pdb -c good.cst -o cr1 >/dev/null 2>&1
check "short form expanded to AtomPair" "$(grep -c '^AtomPair' cr1/input_files/fixed.cst)" 2
check "chain letter kept in the restraint" \
    "$(awk '/^AtomPair/{print $3; exit}' cr1/input_files/fixed.cst)" 2A
# The failure the original could not catch: an atom that is not in the file.
cat > bad.cst <<'EOF'
2A:H    2B:H    1.5  6.5
EOF
"$R/run/rosetta-constrained-relax.sh" -p fixed.pdb -c bad.cst -o cr2 >/dev/null 2>&1
check "missing restrained atom rejected" "$?" 1
printf 'this is not a restraint\n' > junk.cst
"$R/run/rosetta-constrained-relax.sh" -p fixed.pdb -c junk.cst -o cr3 >/dev/null 2>&1
check "unparseable restraint rejected"   "$?" 1

# ---------------------------------------------------------------------------
echo "== interface and subunit checks"
"$R/run/rosetta-symm-dock.sh" -p fixed.pdb -N 3 -o sd >/dev/null 2>&1
check "multi-chain subunit rejected for symmetry" "$?" 1
"$R/run/rosetta-local-dock.sh" -p fixed.pdb -i AB_Z -o ld >/dev/null 2>&1
check "partner naming an absent chain rejected"   "$?" 1
"$R/run/rosetta-local-dock.sh" -p fixed.pdb -i ABC -o ld2 >/dev/null 2>&1
check "partner string without an underscore rejected" "$?" 1

finish
