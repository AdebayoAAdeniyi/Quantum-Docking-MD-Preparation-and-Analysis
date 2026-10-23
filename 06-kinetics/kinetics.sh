#!/usr/bin/env bash
# =============================================================================
#  kinetics.sh - barriers, reaction energies and transition state theory rates
#                from the CSV tables of the QM modules (Gaussian, ORCA, NWChem, xTB)
#
#  usage: kinetics.sh -r reactions.txt -s summary.csv [-t thermo.csv] [-w] [-c 1M] [-o outdir]
#    reactions.txt   one reaction per line:   NAME  reactants  TS  products
#                    with "+" between species, e.g.
#                        h-abstraction   phenol+OH   ts-abs   phenoxyl+H2O
#                        isomerisation   hcn         ts-iso   hnc
#    -s  summary.csv from any extract script (needs columns name, E..., H, G and
#        lowest_freq; the row of each species that has G is used)
#    -t  thermo.csv with several temperatures: rates at every temperature and an
#        Arrhenius fit (Ea, A) over them
#    -w  Wigner tunnelling correction  kappa = 1 + (h c |nu_imag| / kB T)^2 / 24
#    -c  standard state for bimolecular steps: 1atm (as computed) or 1M (default),
#        which adds RT ln(24.46) = 1.89 kcal/mol at 298 K per extra gas-phase
#        species (the usual conversion to 1 mol/L)
#  Eyring equation: k = kappa (kB T / h) exp(-dG#/RT), concentration units (M^1-n s^-1).
#  Outputs: barriers.csv (kJ/mol and kcal/mol), rates.csv, arrhenius.csv, profile plot data.
# =============================================================================
source "$(dirname "$0")/../lib/common.sh"
rx=""; sum=""; thermo=""; wig=0; std=1M; out=kinetics-results
while getopts "r:s:t:wc:o:h" o; do
    case $o in r) rx=$OPTARG ;; s) sum=$OPTARG ;; t) thermo=$OPTARG ;; w) wig=1 ;; c) std=$OPTARG ;; o) out=$OPTARG ;; *) sed -n '2,22p' "$0"; exit 0 ;; esac
done
[[ -f "$rx" && -f "$sum" ]] || die "need -r reactions.txt and -s summary.csv"
mkdir -p "$out"
# species table from the summary: name, E, H, G, lowest frequency, T
awk -F, 'NR==1 { for (i = 1; i <= NF; i++) c[$i] = i
                  e = c["E_scf"] ? c["E_scf"] : (c["E_final"] ? c["E_final"] : c["E_total"])
                  if (!e || !c["G"]) { print "summary table needs an energy column and G" > "/dev/stderr"; exit 1 }
                  next }
         $c["G"] != "" && !(($1) in seen) { seen[$1] = 1
                  print $1, $e, $c["H"], $c["G"], (c["lowest_freq"] ? $c["lowest_freq"] : ""), (c["T_K"] ? $c["T_K"] : 298.15) }' "$sum" > "$out/.species" || exit 1

awk -v wig="$wig" -v std="$std" -v out="$out" '
    BEGIN { KCAL = 627.509474; KJ = 2625.499639; kB = 1.380649e-23; h = 6.62607015e-34; c = 2.99792458e10
            R = 8.314462618; Rk = 1.987204e-3; fb = out "/.b"; fk = out "/.k"; fp = out "/.p" }
    FNR == NR { E[$1] = $2; H[$1] = $3; G[$1] = $4; F[$1] = $5; T[$1] = $6; next }
    /^[[:space:]]*(#|$)/ { next }
    function sum(list, arr,   n, s, i, p) { n = split(list, p, "+"); s = 0
        for (i = 1; i <= n; i++) { if (!(p[i] in G)) { miss = miss " " p[i]; continue } s += arr[p[i]] } return s }
    function count(list,   p) { return split(list, p, "+") }
    {
        name = $1; r = $2; ts = $3; pr = $4; miss = ""
        Er = sum(r, E); Hr = sum(r, H); Gr = sum(r, G); Et = sum(ts, E); Ht = sum(ts, H); Gt = sum(ts, G)
        Ep = sum(pr, E); Hp = sum(pr, H); Gp = sum(pr, G)
        if (miss != "") { printf "warning: %s skipped, no thermochemistry for:%s\n", name, miss > "/dev/stderr"; next }
        TT = T[ts]; nr = count(r); np = count(pr)
        # 1 atm -> 1 M: each species gains RT ln(RT/P) = RT ln 24.46 (298 K); only the change in the number of species matters
        corr = (std == "1M") ? Rk * TT * log(0.082057366 * TT) : 0
        dGa = (Gt - Gr) * KCAL - (nr - 1) * corr; dGr = (Gp - Gr) * KCAL - (nr - np) * corr
        dEa = (Et - Er) * KCAL; dHa = (Ht - Hr) * KCAL; dEr = (Ep - Er) * KCAL; dHr = (Hp - Hr) * KCAL
        nu = F[ts]; kap = 1
        if (wig && nu < 0) { u = h * c * (-nu) / (kB * TT); kap = 1 + u * u / 24 }
        k = kap * kB * TT / h * exp(-dGa * 4184 / (R * TT))
        if (nu >= 0 || nu == "") printf "warning: %s: %s has no imaginary frequency (lowest %s cm-1); is it a transition state?\n", name, ts, nu > "/dev/stderr"
        printf "%s,%s,%s,%s,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.1f,%.2f,%.2f,%.2f,%.2f,%.2f\n", name, r, ts, pr, TT, dEa, dHa, dGa, dEr, dHr, nu, dGr, dEa * 4.184, dHa * 4.184, dGa * 4.184, dGr * 4.184 > fb
        printf "%s,%.2f,%.1f,%.4f,%.4e,%.4e\n", name, TT, nu, kap, k, (k > 0 ? log(2) / k : 0) > fk
        printf "%s 0 0\n%s 1 %.2f\n%s 2 %.2f\n", name, name, dGa, name, dGr > fp
    }' "$out/.species" "$rx"
{ echo "reaction,reactants,TS,products,T_K,dE#_kcal,dH#_kcal,dG#_kcal,dE_r_kcal,dH_r_kcal,nu_imag_cm1,dG_r_kcal,dE#_kJ,dH#_kJ,dG#_kJ,dG_r_kJ"
  cat "$out/.b" 2>/dev/null; } > "$out/barriers.csv"
{ echo "reaction,T_K,nu_imag_cm1,kappa_Wigner,k_TST,half_life_s(first order)"; cat "$out/.k" 2>/dev/null; } > "$out/rates.csv"
awk '{print $1, $2, $3}' "$out/.p" 2>/dev/null > "$out/profile.dat"
rm -f "$out/.b" "$out/.k" "$out/.p"

# rates at several temperatures + Arrhenius fit ln k = ln A - Ea/RT
if [[ -n "$thermo" ]]; then
    [[ -f "$thermo" ]] || die "no $thermo"
    awk -F, -v sumf="$out/.species" 'BEGIN { while ((getline l < sumf) > 0) { split(l, a, " "); E[a[1]] = a[2] } }
        NR==1 { for (i=1;i<=NF;i++) c[$i]=i; next }
        { g = c["G"] ? $c["G"] : (c["G_corr"] && ($1 in E) ? E[$1] + $c["G_corr"] : ""); if (g != "") print $1, $c["T_K"], g }' "$thermo" > "$out/.gT"
    awk -v wig="$wig" -v std="$std" '
        BEGIN { KCAL = 627.509474; kB = 1.380649e-23; h = 6.62607015e-34; R = 8.314462618; Rk = 1.987204e-3; c = 2.99792458e10 }
        FILENAME ~ /species/ { F[$1] = $5; next }
        FILENAME ~ /gT/ { t = sprintf("%.2f", $2); G[$1, t] = $3; temps[t] = 1; next }
        /^[[:space:]]*(#|$)/ { next }
        { n = split($2, r, "+"); m = split($3, ts, "+")
          for (t in temps) { ok = 1; gr = 0; gt = 0
              for (i = 1; i <= n; i++) { if (!((r[i], t) in G)) ok = 0; gr += G[r[i], t] }
              for (i = 1; i <= m; i++) { if (!((ts[i], t) in G)) ok = 0; gt += G[ts[i], t] }
              if (!ok) continue
              corr = (std == "1M") ? Rk * t * log(0.082057366 * t) : 0
              dg = (gt - gr) * KCAL - (n - 1) * corr; kap = 1; nu = F[ts[1]]
              if (wig && nu < 0) { u = h * c * (-nu) / (kB * t); kap = 1 + u * u / 24 }
              printf "%s,%s,%.2f,%.4e\n", $1, t, dg, kap * kB * t / h * exp(-dg * 4184 / (R * t)) } }' \
        "$out/.species" "$out/.gT" "$rx" | sort -t, -k1,1 -k2,2g > "$out/.kt"
    { echo "reaction,T_K,dG#_kcal,k_TST"; cat "$out/.kt"; } > "$out/rates-vs-T.csv"
    # least squares on ln k vs 1/T for each reaction
    awk -F, '{ x = 1 / $2; y = log($4); n[$1]++; sx[$1] += x; sy[$1] += y; sxx[$1] += x * x; sxy[$1] += x * y; syy[$1] += y * y }
        END { print "reaction,n_T,Ea_kJ_mol,A,R2"
              for (k in n) { if (n[k] < 3) continue
                  b = (n[k] * sxy[k] - sx[k] * sy[k]) / (n[k] * sxx[k] - sx[k] ^ 2); a = (sy[k] - b * sx[k]) / n[k]
                  rr = (n[k] * sxy[k] - sx[k] * sy[k]) / sqrt((n[k] * sxx[k] - sx[k] ^ 2) * (n[k] * syy[k] - sy[k] ^ 2))
                  printf "%s,%d,%.2f,%.4e,%.5f\n", k, n[k], -b * 8.314462618 / 1000, exp(a), rr * rr } }' "$out/.kt" > "$out/arrhenius.csv"
    rm -f "$out/.gT" "$out/.kt"
fi
rm -f "$out/.species"
info "kinetics: $out/barriers.csv, $out/rates.csv$( [[ -n "$thermo" ]] && echo ", $out/rates-vs-T.csv, $out/arrhenius.csv")"
