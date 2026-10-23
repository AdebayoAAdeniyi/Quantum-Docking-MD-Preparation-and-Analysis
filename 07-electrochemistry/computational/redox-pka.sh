#!/usr/bin/env bash
# =============================================================================
#  redox-pka.sh - reduction potentials and pKa values from a thermodynamic
#                 (Born-Haber) cycle, for any QM program
#
#  The free energy in solution of each species is
#      G_sol = G_gas + dG_solv,   dG_solv = E_solv - E_gas + dG(1 atm -> 1 M)
#  with G_gas from a gas phase frequency job and E_gas/E_solv single points at the
#  same level and geometry (gas phase and implicit solvent, e.g. SMD).
#
#  usage: redox-pka.sh -k redox|pka -t energies.csv -r reactions.txt [options]
#    energies.csv   name,G_gas,E_gas,E_solv   (Hartree; make it with
#                   gaussian-redox-table.sh or nwchem-redox-table.sh)
#    reactions.txt  redox: NAME  OXIDISED  REDUCED  [n_electrons]
#                   pka:   NAME  ACID      BASE                  (HA -> A- + H+)
#    -R REF      reference electrode: SHE (4.28 V, default), SCE (4.52), Fc (vs Fc+/Fc in
#                acetonitrile, 4.28 + 0.63 = 4.91), or an absolute potential in V
#    -T K        temperature (default 298.15)
#    -s SOLV     water (default) or a value: solvation free energy of H+ in kcal/mol
#    -e          include the free energy of the gas phase electron, -0.867 kcal/mol
#                (Fermi-Dirac statistics, Bartmess, J. Phys. Chem. 1994, 98, 6420)
#  Constants (all can be overridden in the environment):
#    E_abs(SHE) = 4.28 V (Isse and Gennaro, J. Phys. Chem. B 2010, 114, 7894)
#    E(SCE) = +0.241 V vs SHE;  E(Fc+/Fc) = +0.63 V vs SHE in acetonitrile
#      (Pavlishchuk and Addison, Inorg. Chim. Acta 2000, 298, 97)
#    G_gas(H+) = -6.28 kcal/mol, dG_solv(H+, water) = -265.9 kcal/mol, 1 M state
#      (Tissandier et al., J. Phys. Chem. A 1998, 102, 7787; Kelly, Cramer and
#       Truhlar, J. Phys. Chem. B 2006, 110, 16066)
#  Outputs: redox.csv or pka.csv in the current folder (or -o DIR)
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"
kind=""; tab=""; rx=""; ref=SHE; T=298.15; solv=water; elec=0; out=.
while getopts "k:t:r:R:T:s:eo:h" o; do
    case $o in k) kind=$OPTARG ;; t) tab=$OPTARG ;; r) rx=$OPTARG ;; R) ref=$OPTARG ;; T) T=$OPTARG ;;
               s) solv=$OPTARG ;; e) elec=1 ;; o) out=$OPTARG ;; *) sed -n '2,32p' "$0"; exit 0 ;; esac
done
[[ "$kind" =~ ^(redox|pka)$ && -f "$tab" && -f "$rx" ]] || die "need -k redox|pka, -t energies.csv and -r reactions.txt"
mkdir -p "$out"
ESHE=${QDM_E_SHE:-4.28}
case $ref in
    SHE) eref=$ESHE ;; SCE) eref=$(calc "$ESHE + ${QDM_E_SCE:-0.241}") ;;
    Fc|FC|fc) eref=$(calc "$ESHE + ${QDM_E_FC:-0.63}") ;;
    *) [[ "$ref" =~ ^[0-9.]+$ ]] || die "unknown reference $ref"; eref=$ref ;;
esac
[[ $solv == water ]] && gsH=${QDM_DGSOLV_H:--265.9} || gsH=$solv
awk -F, -v kind="$kind" -v T="$T" -v eref="$eref" -v ref="$ref" -v elec="$elec" -v gsH="$gsH" -v GgH="${QDM_GGAS_H:--6.28}" '
    BEGIN { KCAL = 627.509474; F = 23.06054783; R = 1.987204259e-3
            std = R * T * log(0.0820573661 * T)        # 1 atm -> 1 M, kcal/mol (1.89 at 298.15 K)
            ge = elec ? -0.867 : 0 }
    FNR == NR { if (FNR == 1 || $0 ~ /^#/) next
                gsolv[$1] = ($4 - $3) * KCAL + std; ggas[$1] = $2 * KCAL; Gsol[$1] = ggas[$1] + gsolv[$1]; next }
    /^[[:space:]]*(#|$)/ { next }
    { split($0, w, /[ \t]+/); if (w[1] == "") { for (i = 1; i < length(w); i++) w[i] = w[i+1] }
      if (!(w[2] in Gsol) || !(w[3] in Gsol)) { printf "warning: %s: missing %s or %s in the energy table\n", w[1], w[2], w[3] > "/dev/stderr"; next }
      if (kind == "redox") {
          n = (w[4] == "" ? 1 : w[4])
          dGs = Gsol[w[3]] - Gsol[w[2]] - n * ge                 # O + n e- -> R
          dGg = ggas[w[3]] - ggas[w[2]] - n * ge
          E = -dGs / (n * F) - eref
          if (!hdr++) print "name,oxidised,reduced,n,dG_gas_kcal,dG_solv_ox_kcal,dG_solv_red_kcal,dG_red_sol_kcal,dG_red_sol_eV,E_vs_" ref "_V"
          printf "%s,%s,%s,%d,%.2f,%.2f,%.2f,%.2f,%.4f,%.3f\n", w[1], w[2], w[3], n, dGg, gsolv[w[2]], gsolv[w[3]], dGs, dGs / F, E
      } else {
          GH = GgH + gsH + std                                     # G_sol(H+) on the 1 M scale
          dG = Gsol[w[3]] + GH - Gsol[w[2]]
          if (!hdr++) print "name,acid,base,dG_deprot_sol_kcal,pKa,PA_gas_proxy_dG_kcal"
          printf "%s,%s,%s,%.2f,%.2f,%.2f\n", w[1], w[2], w[3], dG, dG / (2.302585093 * R * T), ggas[w[3]] + GgH - ggas[w[2]]
      } }' "$tab" "$rx" > "$out/$kind.csv"
info "$kind: $out/$kind.csv (reference $ref = $eref V absolute)"
