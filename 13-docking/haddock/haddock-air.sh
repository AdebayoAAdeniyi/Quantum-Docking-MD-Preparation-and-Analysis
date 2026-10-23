#!/usr/bin/env bash
# =============================================================================
#  haddock-air.sh - ambiguous interaction restraints (AIR) for HADDOCK
#
#  usage: haddock-air.sh -a "45,47,52" [-p "44,48"] -A "12,15" [-P "11,16"]
#                        [-s A -S B] [-d 2.0] > ambig.tbl
#    -a/-p  active / passive residues of molecule 1 (segid -s, default A)
#    -A/-P  active / passive residues of molecule 2 (segid -S, default B)
#    -d     upper distance bound in A (default 2.0, the HADDOCK default)
#  Every active residue of one partner gets one ambiguous restraint to all
#  active and passive residues of the other, as described in the HADDOCK
#  documentation (Dominguez, Boelens and Bonvin, J. Am. Chem. Soc. 2003, 125, 1731).
#  Active residues usually come from experiment (NMR shifts, mutagenesis) or from
#  predicted interface residues with relative solvent accessibility above 40 %
#  (the NACCESS cut-off used in the original work).
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"
a=""; p=""; A=""; P=""; s=A; S=B; dmax=2.0
while getopts "a:p:A:P:s:S:d:h" o; do case $o in a) a=$OPTARG ;; p) p=$OPTARG ;; A) A=$OPTARG ;; P) P=$OPTARG ;; s) s=$OPTARG ;; S) S=$OPTARG ;; d) dmax=$OPTARG ;; *) sed -n '2,17p' "$0"; exit 0 ;; esac; done
[[ -n "$a" && -n "$A" ]] || die "give the active residues of both partners (-a and -A)"
air() {   # air ACTIVE SEG1 "TARGETS" SEG2
    local r t first
    for r in ${1//,/ }; do
        echo "assign ( resid $r and segid $2 )"
        echo "       ("; first=1
        for t in ${3//,/ }; do ((first)) || echo "     or"; echo "        ( resid $t and segid $4 )"; first=0; done
        echo "       )  $dmax $dmax 0.0"; echo "!"
    done
}
echo "! HADDOCK ambiguous interaction restraints, made by haddock-air.sh"
echo "! molecule $s: active $a${p:+, passive $p}; molecule $S: active $A${P:+, passive $P}"
echo "!"
air "$a" "$s" "$A${P:+,$P}" "$S"
air "$A" "$S" "$a${p:+,$p}" "$s"
