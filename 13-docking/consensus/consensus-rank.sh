#!/usr/bin/env bash
# =============================================================================
#  consensus-rank.sh - combine the scores of several docking programs
#
#  usage: consensus-rank.sh -i scores.csv [-H "glide:low,vina:low,haddock:low,plp:high"]
#                           [-s 10] [-o consensus.csv]
#    scores.csv   ligand,prog1,prog2,...   (header row with program names;
#                 empty, NA or NIL marks a missing score)
#    -H   which direction is better for each program (default: low for all,
#         as for energies; PLP, ChemPLP or GOLD fitness are "high")
#    -s   sigma of the exponential consensus ranking (default 10)
#  Methods (Palacio-Rodriguez, Lans, Cavasotto and Cossio, Sci. Rep. 2019, 9, 5142):
#    AASS  average of auto-scaled scores, (s - best)/(worst - best) -> 0 is best
#    Z     average Z-score, signed so that negative is better
#    RbR   rank by rank: mean rank over programs
#    ECR   exponential consensus ranking: sum over programs of exp(-r/sigma)/sigma
#  Missing scores are skipped for that ligand (the mean is over the programs
#  that docked it), which is what the original scripts did with "NIL".
#  Output is sorted by ECR (highest first), with the rank of each method.
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"
in=""; dirs=""; sigma=10; out=consensus.csv
while getopts "i:H:s:o:h" o; do case $o in i) in=$OPTARG ;; H) dirs=$OPTARG ;; s) sigma=$OPTARG ;; o) out=$OPTARG ;; *) sed -n '2,21p' "$0"; exit 0 ;; esac; done
[[ -f "$in" ]] || die "give the score table with -i"
awk -F, -v dirs="$dirs" -v sig="$sigma" '
    function miss(v) { return (v == "" || v == "NA" || v == "NIL" || v == "nan") }
    NR == 1 { np = NF - 1; for (j = 2; j <= NF; j++) { prog[j] = $j; sgn[j] = 1 }
              n = split(dirs, d, ","); for (i = 1; i <= n; i++) { split(d[i], kv, ":")
                  found = 0; for (j = 2; j <= NF; j++) if (prog[j] == kv[1]) { found = 1; if (kv[2] == "high") sgn[j] = -1 }
                  if (!found) printf "warning: -H names \"%s\", which is not a column of the table\n", kv[1] > "/dev/stderr" }
              next }
    { nl++; lig[nl] = $1; for (j = 2; j <= np + 1; j++) { v[nl, j] = $j; if (!miss($j)) { s = sgn[j] * $j; x[nl, j] = s; c[j]++; sum[j] += s; ss[j] += s * s
          if (!(j in mn) || s < mn[j]) mn[j] = s; if (!(j in mx) || s > mx[j]) mx[j] = s } } }
    END {
        for (j = 2; j <= np + 1; j++) { m[j] = sum[j] / c[j]; sd[j] = sqrt(ss[j] / c[j] - m[j] ^ 2) }
        # rank of every ligand within each program (1 = best = lowest signed score)
        for (j = 2; j <= np + 1; j++) for (i = 1; i <= nl; i++) if ((i, j) in x) { r = 1; for (k = 1; k <= nl; k++) if ((k, j) in x && x[k, j] < x[i, j]) r++; rk[i, j] = r }
        for (i = 1; i <= nl; i++) { a = z = rb = e = 0; k = 0
            for (j = 2; j <= np + 1; j++) if ((i, j) in x) { k++
                a += (mx[j] > mn[j] ? (x[i, j] - mn[j]) / (mx[j] - mn[j]) : 0)
                z += (sd[j] > 0 ? (x[i, j] - m[j]) / sd[j] : 0); rb += rk[i, j]; e += exp(-rk[i, j] / sig) / sig }
            AA[i] = a / k; ZZ[i] = z / k; RB[i] = rb / k; EC[i] = e; NP[i] = k }
        # ranks of the consensus scores themselves
        for (i = 1; i <= nl; i++) { ra = rz = rr = re = 1
            for (k = 1; k <= nl; k++) { if (AA[k] < AA[i]) ra++; if (ZZ[k] < ZZ[i]) rz++; if (RB[k] < RB[i]) rr++; if (EC[k] > EC[i]) re++ }
            line[i] = sprintf("%s,%d,%.4f,%d,%.4f,%d,%.2f,%d,%.5f,%d", lig[i], NP[i], AA[i], ra, ZZ[i], rz, RB[i], rr, EC[i], re)
            for (j = 2; j <= np + 1; j++) line[i] = line[i] "," v[i, j]; key[i] = re }
        printf "ligand,n_programs,AASS,rank_AASS,Zscore,rank_Z,RbR,rank_RbR,ECR,rank_ECR"; for (j = 2; j <= np + 1; j++) printf ",%s", prog[j]; print ""
        for (r = 1; r <= nl; r++) for (i = 1; i <= nl; i++) if (key[i] == r) print line[i]
    }' "$in" > "$out"
info "consensus ranking: $out"
