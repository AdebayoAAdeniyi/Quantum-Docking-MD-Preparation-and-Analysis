#!/usr/bin/env bash
# =============================================================================
#  docking-best-poses.sh - best pose per ligand from multi-pose SDF output of
#                          Glide, MOE, GOLD or any program that writes scores
#                          as SDF data fields
#
#  usage: docking-best-poses.sh -i poses.sdf -f FIELD [-n NAMEFIELD] [-H] [-o out]
#    -f   score field: r_i_docking_score (Glide), S (MOE), Gold.PLP.Fitness ...
#    -n   field holding the ligand name (default: the title line of each record;
#         MOE writes "mol" or "molecule name", Glide keeps the title)
#    -H   higher score is better (GOLD fitness); default lower is better
#  Writes OUT-best.sdf (one record per ligand), OUT-best.csv (sorted) and
#  OUT-all.csv (every pose). Replaces the MOE mol2/sdf splitting scripts and the
#  Glide log parsing of the original work; the files are read with awk only, so
#  no Schrodinger or MOE licence is needed for the analysis.
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"
in=""; fld=""; nf=""; high=0; out=""
while getopts "i:f:n:Ho:h" o; do case $o in i) in=$OPTARG ;; f) fld=$OPTARG ;; n) nf=$OPTARG ;; H) high=1 ;; o) out=$OPTARG ;; *) sed -n '2,17p' "$0"; exit 0 ;; esac; done
[[ -f "$in" && -n "$fld" ]] || die "need -i poses.sdf and -f score field"
out=${out:-$(strip_ext "$in")}
awk -v fld="$fld" -v nf="$nf" -v high="$high" -v out="$out" '
    BEGIN { RS = "\\$\\$\\$\\$\r?\n"; ORS = "" }
    NF == 0 { next }
    { rec = $0; n = split(rec, L, "\n"); title = L[1]; nm = title; sc = ""
      for (i = 1; i < n; i++) {
          if (L[i] ~ "^> *<" fld ">") sc = L[i+1]
          if (nf != "" && L[i] ~ "^> *<" nf ">") nm = L[i+1] }
      gsub(/\r/, "", sc); gsub(/\r/, "", nm); gsub(/^ +| +$/, "", nm)
      if (sc == "") next
      np++; print nm "," sc "\n" > (out "-all.csv")
      if (!(nm in best) || (high ? sc + 0 > best[nm] + 0 : sc + 0 < best[nm] + 0)) { best[nm] = sc; brec[nm] = rec } }
    END { for (k in best) { print brec[k] "$$$$\n" > (out "-best.sdf"); print k "," best[k] "\n" > (out ".tmp") } }' "$in"
[[ -f "$out.tmp" ]] || die "no records with field <$fld> in $in"
{ echo "ligand,$fld"; if ((high)); then sort -t, -k2,2gr "$out.tmp"; else sort -t, -k2,2g "$out.tmp"; fi; } > "$out-best.csv"
sed -i "1i ligand,$fld" "$out-all.csv"; rm -f "$out.tmp"
info "$(($(wc -l < "$out-best.csv") - 1)) ligands, $(($(wc -l < "$out-all.csv") - 1)) poses -> $out-best.csv, $out-best.sdf"
