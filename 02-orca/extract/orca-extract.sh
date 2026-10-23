#!/usr/bin/env bash
# =============================================================================
#  orca-extract.sh - collect ORCA results of many molecules into CSV tables
#
#  usage: orca-extract.sh -k KIND [-l list] [-p parent] [-o outdir] [-r refs | -R ref.out]
#   KIND      table
#   summary   status, optimisation converged, energies, ZPE, H, G, HOMO/LUMO,
#             dipole, imaginary frequencies (one row per step output)
#   thermo    one row per temperature (temps step)
#   charges   Mulliken, Loewdin, Hirshfeld, CHELPG charges, Mulliken spin
#   freq      frequencies, IR intensities, Raman activities
#   uv        TD-DFT states with main orbital transitions
#   socuv     spin-orbit corrected spectrum        socme  |<T|HSO|S>| couplings
#   casscf    CASSCF and NEVPT2 transition energies, active space occupations
#   neb       energy profile of the NEB path (CI and TS marked)
#   esd       fluorescence / phosphorescence rates with FC and HT shares
#   nmr       shieldings (+ shifts with -r "C=..,H=.." or -R reference.out)
#   geometry  final structures as NAME-STEP-final.xyz
#   all       everything above
#  Outputs are found as PARENT/NAME/NAME-*.out (default parent orca-jobs).
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"
AWKO="$QDM_ROOT/lib/awk/orca.awk"
kind=""; list=""; parent=orca-jobs; out=orca-results; refs=""; reflog=""
while getopts "k:l:p:o:r:R:h" o; do
    case $o in k) kind=$OPTARG ;; l) list=$OPTARG ;; p) parent=$OPTARG ;; o) out=$OPTARG ;;
               r) refs=$OPTARG ;; R) reflog=$OPTARG ;; *) sed -n '2,21p' "$0"; exit 0 ;; esac
done
shift $((OPTIND-1))
[[ -n "$kind" ]] || die "choose a table with -k (see -h)"
if [[ -n "$list" ]]; then mapfile -t names < <(read_list "$list")
else mapfile -t names < <(find "$parent" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort); fi
outs=()
for n in "${names[@]}"; do
    found=0
    for f in "$parent/$n/$n"-*.out "$parent/$n/$n.out"; do [[ -f "$f" ]] && { outs+=("$f"); found=1; }; done
    ((found)) || warn "$n: no ORCA output in $parent/$n"
done
outs+=("$@")
((${#outs[@]})) || die "no ORCA outputs found"
mkdir -p "$out"

table() {
    local k=$1 csv="$out/$1.csv" first=1 f name step
    : > "$csv"
    for f in "${outs[@]}"; do
        name=$(basename "$(dirname "$f")"); step=$(strip_ext "$f"); step=${step#"$name"}; step=${step#-}; step=${step:-main}
        awk -v mode="$k" -v name="$name,$step" -v header=$first -f "$AWKO" "$f" |
            { if ((first)); then sed '1s/^name,/name,step,/'; else cat; fi; } >> "$csv"
        [[ -s "$csv" ]] && first=0
    done
    if [[ $k != summary && $(wc -l < "$csv") -le 1 ]]; then warn "no $k data found"; rm -f "$csv"; return; fi
    info "$k: $csv ($(($(wc -l < "$csv") - 1)) rows)"
}
nmr_shifts() {
    local csv="$out/nmr.csv" spec="$refs"; [[ -f "$csv" ]] || return 0
    if [[ -n "$reflog" ]]; then
        spec=$(awk -v mode=nmr -f "$AWKO" "$reflog" | awk -F, '{s[$3]+=$4; n[$3]++} END{for (e in s) printf "%s=%.4f,", e, s[e]/n[e]}')
        spec=${spec%,}; info "NMR references from $(basename "$reflog"): $spec"
    fi
    [[ -z "$spec" ]] && return 0
    awk -F, -v OFS=, -v spec="$spec" '
        BEGIN { n = split(spec, a, ","); for (i = 1; i <= n; i++) { split(a[i], kv, "="); split(kv[2], ab, ":"); A[kv[1]] = ab[1]; B[kv[1]] = (ab[2] == "" ? 1 : ab[2]) } }
        NR == 1 { print $0, "shift_ppm"; next }
        { e = $4; print $0, ((e in A) ? sprintf("%.2f", A[e] - B[e] * $5) : "") }' "$csv" > "$csv.tmp" && mv "$csv.tmp" "$csv"
}
geometry() {
    local f name step
    for f in "${outs[@]}"; do
        grep -q "CARTESIAN COORDINATES (ANGSTROEM)" "$f" || continue
        name=$(basename "$(dirname "$f")"); step=$(strip_ext "$f")
        awk -v mode=geometry -v name="$step" -f "$AWKO" "$f" > "$out/$step-final.xyz"
    done
    info "geometry: $out/*-final.xyz"
}
case $kind in
    all) for k in summary thermo charges freq uv socuv socme casscf neb esd nmr; do table $k; done; nmr_shifts; geometry ;;
    nmr) table nmr; nmr_shifts ;;
    geometry) geometry ;;
    summary|thermo|charges|freq|uv|socuv|socme|casscf|neb|esd) table "$kind" ;;
    *) die "unknown kind: $kind" ;;
esac
