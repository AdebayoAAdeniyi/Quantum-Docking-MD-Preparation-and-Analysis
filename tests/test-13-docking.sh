#!/usr/bin/env bash
# Docking: a real Vina redocking of the 1IEP crystal ligand, plus the consensus
# ranking and the pose/restraint helpers on their file formats.
source "$(dirname "$0")/testlib.sh"
D=$QDM_ROOT/13-docking
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT; cd "$W"
echo "== consensus ranking"
cat > sc.csv <<'EOT'
ligand,vina,glide,plp
A,-13.2,-11.4,78.1
B,-12.0,-10.9,74.3
C,-14.2,-12.1,81.6
D,-8.6,NIL,41.2
EOT
"$D/consensus/consensus-rank.sh" -i sc.csv -H "plp:high" -o cons.csv 2>/dev/null
check "best ligand first"      "$(awk -F, 'NR==2{print $1}' cons.csv)" C
check "missing score skipped"  "$(awk -F, '$1=="D"{print $2}' cons.csv)" 2
check "AASS of the best = 0"   "$(awk -F, 'NR==2{print $3}' cons.csv)" 0.0000
check "ECR ranks agree"        "$(awk -F, 'NR>1 && $4!=$10{n++} END{print n+0}' cons.csv)" 0
# a fitness score alone: declared "high" the best ligand is C, otherwise it inverts to D
awk -F, 'NR==1{print "ligand,plp"; next} {print $1 "," $4}' sc.csv > plp.csv
"$D/consensus/consensus-rank.sh" -i plp.csv -H "plp:high" -o p1.csv 2>/dev/null
"$D/consensus/consensus-rank.sh" -i plp.csv -o p2.csv 2>/dev/null
check "high score is best"     "$(awk -F, 'NR==2{print $1}' p1.csv)" C
check "direction matters"      "$(awk -F, 'NR==2{print $1}' p2.csv)" D
echo "== HADDOCK restraints"
"$D/haddock/haddock-air.sh" -a 45,47 -p 44 -A 12 -P 11,16 > air.tbl 2>/dev/null
check "one block per active residue" "$(grep -c '^assign' air.tbl)" 3
check "passive residues included"    "$(grep -c 'resid 16 and segid B' air.tbl)" 2
check "distance bounds"              "$(grep -c '2.0 2.0 0.0' air.tbl)" 3
echo "== best pose from a multi-pose SDF"
for i in 1 2 3; do printf 'LIG%d\n\n\n  0  0\nM  END\n> <r_i_docking_score>\n-%d.5\n\n$$$$\n' $((i%2+1)) $i; done > poses.sdf
"$D/glide-moe/docking-best-poses.sh" -i poses.sdf -f r_i_docking_score -o g 2>/dev/null
check "two ligands, best kept" "$(($(wc -l < g-best.csv) - 1))" 2
check "lowest score wins"      "$(awk -F, '$1=="LIG2"{print $2}' g-best.csv)" -3.5
echo "== Vina redocking of 1IEP (needs vina, meeko, obabel)"
need vina; need obabel
command -v mk_prepare_ligand.py >/dev/null || { echo "  Meeko not installed, docking skipped"; finish; exit; }
fetch "https://raw.githubusercontent.com/ccsb-scripps/AutoDock-Vina/develop/example/basic_docking/data/1iep_receptorH.pdb" 1iep_receptorH.pdb || exit 0
fetch "https://raw.githubusercontent.com/ccsb-scripps/AutoDock-Vina/develop/example/basic_docking/data/1iep_ligand.sdf" 1iep_ligand.sdf || exit 0
"$D/vina/vina-prepare.sh" -r "$DATA/1iep_receptorH.pdb" -x "$DATA/1iep_ligand.sdf" -p 6 "$DATA/1iep_ligand.sdf" 2>/dev/null
check "receptor pdbqt"   "$([[ -s vina-setup/receptor.pdbqt ]] && echo yes)" yes
check "box from ligand"  "$(awk -F= '/center_x/{printf "%.1f", $2}' vina-setup/box.txt)" 15.1 0.2
"$D/vina/vina-run.sh" -e 32 -n 9 -c "${CORES:-2}" 2>/dev/null
"$D/vina/vina-extract.sh" -x "$DATA/1iep_ligand.sdf" 2>/dev/null
check "affinity is reasonable" "$(awk -F, 'NR==2{print $2}' vina-results/vina-best.csv)" -12.5 2.5
rms=$(awk -F, 'NR==2{gsub(/"/,"",$3); print $3}' vina-results/vina-best.csv)
check "redocks within 2 A"     "$(awk -v r="$rms" 'BEGIN{print (r + 0 > 0 && r + 0 < 2) ? "yes" : "no"}')" yes
finish
