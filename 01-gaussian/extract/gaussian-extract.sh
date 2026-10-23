#!/usr/bin/env bash
# =============================================================================
#  gaussian-extract.sh - collect results of many Gaussian jobs into CSV tables
#
#  Every log is split into its --Link1-- steps first, so gas phase and solvent
#  energies, optimisation and property steps are reported separately (the old
#  scripts took "the last match" and mixed them). Each row carries the step.
#
#  usage: gaussian-extract.sh -k KIND [options] [-l list -p parent | logs...]
#   KIND       table (one CSV per kind, written to -o DIR)
#   summary    status, charge, mult, energies, ZPE, H, G, Cv, S, HOMO, LUMO,
#              gap, dipole, imaginary frequencies (one row per step)
#   thermo     thermochemistry, one row per temperature (gaussian-followup -k temps)
#   charges    Mulliken, APT, NPA, Hirshfeld, CM5 charges and spin per atom
#   uv         excited states: eV, nm, f, main transitions (% = 2c^2 x 100)
#   freq       frequencies with IR intensities and Raman activities
#   nmr        isotropic shielding, and chemical shifts when -r or -R is given
#   polar      alpha (iso, aniso) and beta (components, total, esu)
#   hyperfine  isotropic Fermi contact couplings (open shell)
#   jcoupling  NMR spin-spin couplings J with FC, SD, PSO, DSO terms
#              (all pairs with |J| >= 0.5 Hz, or only the pairs given by -b)
#   fukui      condensed Fukui functions + IP, EA, mu, eta, omega
#              (from NAME-fukui.log made by gaussian-input.sh -F)
#   geometry   final geometry of each molecule as NAME-final.xyz
#   all        every table above
#  options:
#    -l LIST -p DIR   molecules in DIR/NAME/ (default gaussian-jobs)
#    -o DIR           output folder (default gaussian-results)
#    -r REFS          NMR references, "C=191.12,H=31.42"  (delta = ref - sigma)
#                     or scaled "C=175.7:0.963"           (delta = a - b*sigma)
#    -R TMS.log       take reference shieldings from a reference calculation
#                     done at the same level of theory (recommended)
#    -b PAIRS         atom pairs for jcoupling, e.g. "1-2,3-7"
#  example:
#    gaussian-extract.sh -k all -l list -R tms/tms.log
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"
AWKG="$QDM_ROOT/lib/awk/gaussian.awk"
SPLIT="$(dirname "$0")/split-link1-log.sh"
kind=""; list=""; parent=gaussian-jobs; out=gaussian-results; refs=""; reflog=""; pairs=""
while getopts "k:l:p:o:r:R:b:h" o; do
    case $o in k) kind=$OPTARG ;; l) list=$OPTARG ;; p) parent=$OPTARG ;; o) out=$OPTARG ;;
               r) refs=$OPTARG ;; R) reflog=$OPTARG ;; b) pairs=$OPTARG ;; *) sed -n "2,36p" "$0"; exit 0 ;; esac
done
shift $((OPTIND-1))
[[ -n "$kind" ]] || die "choose a table with -k (see -h)"
# ---- collect logs: NAME -> every NAME*.log in parent/NAME ----
logs=("$@")
if [[ -n "$list" ]]; then
    while read -r n; do
        found=0
        for f in "$parent/$n/$n.log" "$parent/$n/$n"-*.log; do
            [[ -f "$f" ]] && { logs+=("$f"); found=1; }
        done
        ((found)) || warn "$n: no log in $parent/$n"
    done < <(read_list "$list")
fi
((${#logs[@]})) || die "no log files (use -l list -p parent or give files)"
mkdir -p "$out"
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT

# steps LOG -> step files (in a temp folder)
steps() { "$SPLIT" "$1" "$tmp/$(strip_ext "$1")" 2>/dev/null; }

# table KIND : one CSV over all logs and steps, columns "name,step,..."
table() {
    local k=$1 csv="$out/$1.csv" first=1 log st name label
    : > "$csv"
    for log in "${logs[@]}"; do
        [[ "$(basename "$log")" == *-fukui.log ]] && continue      # handled by kind fukui
        name=$(strip_ext "$log")
        while read -r st; do
            label=$(basename "$st" .log); label=${label#"$name".}
            awk -v mode="$k" -v name="$name,$label" -v header=$first -v pairs="$pairs" -f "$AWKG" "$st" |
                { if ((first)); then sed '1s/^name,/name,step,/'; else cat; fi; } >> "$csv"
            [[ -s "$csv" ]] && first=0
        done < <(steps "$log")
    done
    # drop data-less rows of tables that only make sense for some steps
    case $k in
        thermo|uv|freq|nmr|hyperfine|charges|jcoupling)
            [[ $(wc -l < "$csv") -le 1 ]] && { warn "no $k data found"; rm -f "$csv"; return; } ;;
        polar) awk -F, 'NR==1 || $3 != ""' "$csv" > "$tmp/p" && mv "$tmp/p" "$csv" ;;
    esac
    info "$k: $csv ($(($(wc -l < "$csv") - 1)) rows)"
}

# NMR shifts: delta = ref - sigma, or a - b*sigma for scaled references
nmr_shifts() {
    local csv="$out/nmr.csv"; [[ -f "$csv" ]] || return 0
    local spec="$refs"
    if [[ -n "$reflog" ]]; then
        [[ -f "$reflog" ]] || die "reference log not found: $reflog"
        # mean shielding of each element in the reference molecule (e.g. TMS)
        spec=$(awk -v mode=nmr -f "$AWKG" "$reflog" |
               awk -F, '{s[$3]+=$4; n[$3]++} END{for (e in s) printf "%s=%.4f,", e, s[e]/n[e]}')
        spec=${spec%,}
        info "NMR references from $(basename "$reflog"): $spec"
    fi
    [[ -z "$spec" ]] && { info "nmr: shieldings only (give -r or -R for chemical shifts)"; return 0; }
    awk -F, -v OFS=, -v spec="$spec" '
        BEGIN { n = split(spec, a, ",")
                for (i = 1; i <= n; i++) { split(a[i], kv, "="); split(kv[2], ab, ":")
                    A[kv[1]] = ab[1]; B[kv[1]] = (ab[2] == "" ? 1 : ab[2]) } }
        NR == 1 { print $0, "shift_ppm", "reference"; next }
        { e = $4; if (e in A) printf "%s,%.2f,%s\n", $0, A[e] - B[e] * $5, A[e] (B[e] != 1 ? ":" B[e] : "")
          else print $0, "", "" }' "$csv" > "$tmp/n" && mv "$tmp/n" "$csv"
}

# Conceptual DFT from the three steps of NAME-fukui.log (N, N+1, N-1)
fukui() {
    local csvA="$out/fukui-atoms.csv" csvG="$out/fukui-global.csv" log name s n first=1
    echo "name,IP_eV,EA_eV,mu_eV,eta_eV,softness_1_eV,omega_eV,E_N,E_N+1,E_N-1" > "$csvG"
    : > "$csvA"
    for log in "${logs[@]}"; do
        [[ "$(basename "$log")" == *-fukui.log ]] || continue
        name=$(basename "$log" -fukui.log)
        mapfile -t s < <(steps "$log")
        ((${#s[@]} == 3)) || { warn "$name: expected 3 steps (N, N+1, N-1) in $log, found ${#s[@]}"; continue; }
        for n in 0 1 2; do
            grep -q "Normal termination" "${s[$n]}" || { warn "$name: step $((n+1)) of the Fukui job did not finish"; continue 2; }
            awk -v mode=charges -f "$AWKG" "${s[$n]}" > "$tmp/q$n"
            awk -v mode=summary -f "$AWKG" "${s[$n]}" | cut -d, -f6 > "$tmp/e$n"
        done
        # global descriptors (Parr and Pearson 1983; Parr, Szentpaly and Liu 1999)
        paste -d, "$tmp/e0" "$tmp/e1" "$tmp/e2" | awk -F, -v n="$name" -v h="$HARTREE2EV" '{
            ip = ($3 - $1) * h; ea = ($1 - $2) * h; mu = -(ip + ea) / 2; eta = (ip - ea) / 2
            printf "%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%s,%s,%s\n", n, ip, ea, mu, eta, 1/(2*eta), mu*mu/(2*eta), $1, $2, $3 }' >> "$csvG"
        # condensed Fukui functions from NPA and Hirshfeld charges
        paste -d, "$tmp/q0" "$tmp/q1" "$tmp/q2" | awk -F, -v n="$name" -v hdr=$first '
            BEGIN { if (hdr) print "name,atom,element,q_N_npa,f+_npa,f-_npa,f0_npa,q_N_hirsh,f+_hirsh,f-_hirsh,f0_hirsh,dual_descriptor_hirsh" }
            # columns per block: name,atom,el,mul,apt,npa,hir,cm5,spin -> npa=6, hir=7 ; blocks start at 1, 10, 19
            { p0=$6; p1=$15; p2=$24; h0=$7; h1=$16; h2=$25
              fp=p0-p1; fm=p2-p0; hp=h0-h1; hm=h2-h0
              npa = (p0 == "" || p1 == "" || p2 == "") ? ",,," : sprintf("%.4f,%.4f,%.4f", fp, fm, (p2-p1)/2)
              hir = (h0 == "" || h1 == "" || h2 == "") ? ",,,," : sprintf("%.4f,%.4f,%.4f,%.4f", hp, hm, (h2-h1)/2, hp-hm)
              printf "%s,%s,%s,%s,%s,%s,%s\n", n, $2, $3, p0, npa, h0, hir }' >> "$csvA"
        first=0
    done
    info "fukui: $csvG and $csvA"
}

geometry() {
    local log name st
    for log in "${logs[@]}"; do
        [[ "$(basename "$log")" == *-* && "$(basename "$log")" != "$(basename "$(dirname "$log")").log" ]] && continue
        name=$(strip_ext "$log")
        st=$(steps "$log" | grep -E -- '-(opt|ts|qst)' | tail -n 1)
        [[ -z "$st" ]] && st=$(steps "$log" | tail -n 1)
        awk -v mode=geometry -v name="$name" -f "$AWKG" "$st" > "$out/$name-final.xyz"
    done
    info "geometry: $out/*-final.xyz"
}

case $kind in
    all) for k in summary thermo charges uv freq nmr jcoupling polar hyperfine; do table $k; done
         nmr_shifts; fukui; geometry ;;
    nmr) table nmr; nmr_shifts ;;
    fukui) fukui ;;
    geometry) geometry ;;
    summary|thermo|charges|uv|freq|polar|hyperfine|jcoupling) table "$kind" ;;
    *) die "unknown kind: $kind" ;;
esac
