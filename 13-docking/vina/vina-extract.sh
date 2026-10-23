#!/usr/bin/env bash
# =============================================================================
#  vina-extract.sh - scores, best poses and redocking RMSD
#  usage: vina-extract.sh [-s vina-setup] [-x reference_ligand] [-o vina-results]
#  Writes vina-scores.csv (all poses), vina-best.csv (pose 1 per ligand, sorted),
#  best poses as SDF (mk_export.py) and, with -x, the heavy atom RMSD of every
#  ligand's pose 1 against the reference (symmetry-aware, Open Babel obrms).
#  A redocked crystal ligand within 2 A of its crystal pose is the usual check
#  that the box and settings reproduce the experimental binding mode.
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"
s=vina-setup; ref=""; out=vina-results
while getopts "s:x:o:h" o; do case $o in s) s=$OPTARG ;; x) ref=$OPTARG ;; o) out=$OPTARG ;; *) sed -n '2,10p' "$0"; exit 0 ;; esac; done
mkdir -p "$out/poses"
echo "ligand,pose,affinity_kcal_mol,rmsd_lb,rmsd_ub" > "$out/vina-scores.csv"
for f in "$s"/docked/*_out.pdbqt; do
    [[ -f "$f" ]] || continue
    n=$(basename "$f" _out.pdbqt)
    awk -v n="$n" '/^REMARK VINA RESULT:/ { p++; printf "%s,%d,%s,%s,%s\n", n, p, $4, $5, $6 }' "$f" >> "$out/vina-scores.csv"
    if command -v mk_export.py >/dev/null; then mk_export.py "$f" -s "$out/poses/$n.sdf" >/dev/null 2>&1 || obabel "$f" -O "$out/poses/$n.sdf" 2>/dev/null
    else obabel "$f" -O "$out/poses/$n.sdf" 2>/dev/null; fi
done
{ echo "ligand,best_affinity_kcal_mol$( [[ -n "$ref" ]] && echo ',rmsd_to_reference_A')"
  awk -F, 'NR>1 && $2==1 {print $1 "," $3}' "$out/vina-scores.csv" | sort -t, -k2,2g | while IFS=, read -r n e; do
      if [[ -n "$ref" && -s "$out/poses/$n.sdf" ]]; then
          obabel "$out/poses/$n.sdf" -l 1 -O "$out/poses/$n-top.sdf" 2>/dev/null
          r=$(obrms -f "$ref" "$out/poses/$n-top.sdf" 2>/dev/null | awk '{print $NF}' | head -1)
          # obrms returns inf when the two molecules are not the same compound
          [[ -z "$r" || "$r" == inf || "$r" == nan ]] && r="not the reference compound"
          echo "$n,$e,\"$r\""
      else echo "$n,$e"; fi
  done; } > "$out/vina-best.csv"
info "Vina: $out/vina-best.csv ($(($(wc -l < "$out/vina-best.csv") - 1)) ligands)"
cat "$out/vina-best.csv" >&2
