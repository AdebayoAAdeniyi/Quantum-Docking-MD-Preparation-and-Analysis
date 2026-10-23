#!/usr/bin/env bash
# =============================================================================
#  orca-esd-input.sh - fluorescence and phosphorescence rate jobs (ORCA ESD module)
#
#  Run after the "excited" step of orca-input.sh has finished, i.e. when every
#  folder holds NAME-GS.out/.hess/.xyz (S0) and NAME-TG.out/.hess (T1).
#
#   * fluorescence:   ESD(FLUOR) from S1, with the excited state Hessian
#                     approximated from the ground state one (HESSFLAG AHAS)
#   * phosphorescence: ESD(PHOSP) for the three T1 sublevels (IROOT 1, 2, 3 in
#                     three $new_job steps), with the S0 and T1 Hessians and the
#                     adiabatic S0-T1 energy DELE, which ORCA requires as input
#  DELE = E(T1) - E(S0) at their own optimised geometries, electronic energies
#  only (add -Z to include the zero point energies), in cm-1.
#
#  usage: orca-esd-input.sh [-l list] [-p parent] [-s smd:solvent] [-Z]
#                           [-m METHOD] [-b BASIS] [-t nroots] [-w linewidth]
#  Use the same method and basis as the excited step (defaults B3LYP, def2-TZVP).
#  The original workflow used THF for luminescence of biological probes, since
#  emission is usually measured in a hydrophobic environment.
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"
source "$QDM_ROOT/lib/jobscript.sh"
load_cluster_conf
AWKO="$QDM_ROOT/lib/awk/orca.awk"
list=""; parent=orca-jobs; solv=""; zpe=0; method=B3LYP; basis=def2-TZVP; nroots=25; linew=75
nproc=${CORES:-8}; mem=${MEM:-16G}
while getopts "l:p:s:Zm:b:t:w:n:M:h" o; do
    case $o in l) list=$OPTARG ;; p) parent=$OPTARG ;; s) solv=$OPTARG ;; Z) zpe=1 ;; m) method=$OPTARG ;;
               b) basis=$OPTARG ;; t) nroots=$OPTARG ;; w) linew=$OPTARG ;; n) nproc=$OPTARG ;; M) mem=$OPTARG ;;
               *) sed -n '2,22p' "$0"; exit 0 ;; esac
done
CORES=$nproc; MEM=$mem
memmb=$(echo "$mem" | awk '{v=$1+0; u=toupper($1); if (u ~ /T/) v*=1048576; else if (u ~ /G/) v*=1024; print v}')
maxcore=$(( memmb * 3 / 4 / nproc ))
solkw=""; solblock=""
if [[ -n "$solv" ]]; then
    sname=${solv#*:}; solkw=" CPCM($sname)"
    [[ "${solv%%:*}" == smd ]] && solblock=$'%cpcm\n  smd true\n  SMDsolvent "'"$sname"$'"\nend'
fi
if [[ -n "$list" ]]; then mapfile -t names < <(read_list "$list")
else mapfile -t names < <(find "$parent" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort); fi

energy() {   # total energy (+ ZPE with -Z) of an ORCA output, in Eh
    awk -v mode=summary -f "$AWKO" "$1" | awk -F, -v z=$zpe '{ printf "%.10f", $7 + (z ? $8 : 0) }'
}
head_esd() {
    echo "! $method $basis TightSCF$solkw $1"
    echo "%pal nprocs $nproc end"; echo "%maxcore $maxcore"
    [[ -n "$solblock" ]] && echo "$solblock"
    echo "%scf MaxIter 500 end"
    return 0
}
for n in "${names[@]}"; do
    d="$parent/$n"
    for f in "$n-GS.out" "$n-TG.out" "$n-GS.hess" "$n-TG.hess" "$n-GS.xyz"; do
        [[ -f "$d/$f" ]] || { warn "$n: $f missing, skipped"; continue 2; }
    done
    grep -q "ORCA TERMINATED NORMALLY" "$d/$n-GS.out" "$d/$n-TG.out" || { warn "$n: S0 or T1 job did not finish"; continue; }
    q=$(awk -v mode=summary -f "$AWKO" "$d/$n-GS.out" | cut -d, -f4)
    egs=$(energy "$d/$n-GS.out"); etg=$(energy "$d/$n-TG.out")
    dele=$(awk -v a="$etg" -v b="$egs" 'BEGIN{printf "%.1f", (a - b) * 219474.6313632}')
    awk -v x="$dele" 'BEGIN{exit !(x <= 0)}' && { warn "$n: T1 below S0 (DELE = $dele cm-1), check the jobs"; continue; }
    cp "$d/$n-GS.xyz" "$d/$n-esd.xyz"
    {
        head_esd "ESD(FLUOR)"
        printf '%%tddft\n  NRoots 5\n  IRoot 1\nend\n'
        printf '%%esd\n  GSHessian "%s"\n  HessFlag AHAS\n  DoHT true\n  Lines Voigt\n  LineW %s\n  InLineW 200\nend\n' "$n-GS.hess" "$linew"
        echo "* xyzfile $q 1 $n-esd.xyz"
    } > "$d/$n-fluor.inp"
    {
        for r in 1 2 3; do
            ((r > 1)) && echo '$new_job'
            head_esd "ESD(PHOSP) RI-SOMF(1X)"
            printf '%%tddft\n  NRoots %s\n  DoSOC true\n  TDA false\n  IRoot %s\nend\n' "$nroots" "$r"
            printf '%%esd\n  GSHessian "%s"\n  TSHessian "%s"\n  DoHT true\n  LineW 200\n  DELE %s\nend\n' "$n-GS.hess" "$n-TG.hess" "$dele"
            echo "* xyzfile $q 1 $n-esd.xyz"
        done
    } > "$d/$n-phosp.inp"
    write_jobscript "esd-$n" "$PWD/$d" ORCA_SETUP "ORCA_EXE=\"${ORCA_EXE:-orca}\"" \
        "\$ORCA_EXE $n-fluor.inp > $n-fluor.out" "\$ORCA_EXE $n-phosp.inp > $n-phosp.out" > "$d/run-esd.sh"
    chmod +x "$d/run-esd.sh"
    info "$n: DELE = $dele cm-1 ($(awk -v x="$dele" 'BEGIN{printf "%.3f", x/8065.544}') eV); submit with qdm-submit -s run-esd.sh $d"
done
