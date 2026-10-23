#!/usr/bin/env bash
# =============================================================================
#  plumed-input.sh - PLUMED input for well-tempered metadynamics, optionally
#                    inside a funnel restraint (funnel metadynamics)
#
#  usage: plumed-input.sh -p complex.pdb -l LIGAND_SEL [-a "SEL_A" -b "SEL_B"] [options] > plumed.dat
#    Selections are residue numbers ("315,381"), residue ranges ("40-45") or a
#    residue name (LIG); atoms are taken from the PDB (same numbering as the MD).
#    -l  ligand
#    -a  pocket atoms that define point A, deep in the binding site
#    -b  atoms that define point B, towards the solvent (the funnel axis is A -> B)
#    CVs (default): projection of the ligand centre on the A->B axis (fps.lp)
#    -c TYPE  cv: "funnel" (default) or "distance" (ligand - pocket centre distance)
#    -z       zcc, where the cone turns into the cylinder, nm (default 1.8)
#    -r       cylinder radius, nm (default 0.1)      -A  cone half-angle, rad (default 0.6)
#    -m/-M    lower/upper walls on the projection, nm (default -0.2 / 3.0)
#    -s SIGMA -H HEIGHT (kJ/mol) -g BIASFACTOR -T TEMP -P PACE (steps)
#  The funnel is built from standard PLUMED actions (PROJECTION_ON_AXIS, CUSTOM,
#  UPPER_WALLS), so the optional funnel module does not have to be compiled in.
#  Method: Limongelli, Bonomi and Parrinello, PNAS 2013, 110, 6358.
#  Amber reads the file with  plumed=1, plumedfile='plumed.dat'  in &cntrl;
#  GROMACS needs a PLUMED-patched build (mdrun -plumed plumed.dat); CP2K: cp2k-input.sh -P.
# =============================================================================
source "$(dirname "$0")/../lib/common.sh"
pdb=""; lig=""; sa=""; sb=""; cv=funnel; zcc=1.8; rcyl=0.1; alpha=0.6; zmin=-0.2; zmax=3.0
sigma=0.05; height=1.2; gamma=10; T=300; pace=500
while getopts "p:l:a:b:c:z:r:A:m:M:s:H:g:T:P:h" o; do
    case $o in p) pdb=$OPTARG ;; l) lig=$OPTARG ;; a) sa=$OPTARG ;; b) sb=$OPTARG ;; c) cv=$OPTARG ;; z) zcc=$OPTARG ;; r) rcyl=$OPTARG ;;
               A) alpha=$OPTARG ;; m) zmin=$OPTARG ;; M) zmax=$OPTARG ;; s) sigma=$OPTARG ;; H) height=$OPTARG ;; g) gamma=$OPTARG ;;
               T) T=$OPTARG ;; P) pace=$OPTARG ;; *) sed -n '2,25p' "$0" >&2; exit 0 ;; esac
done
[[ -f "$pdb" && -n "$lig" && -n "$sa" ]] || die "need -p PDB, -l ligand and -a pocket selection"
atoms() {   # atoms SELECTION -> comma separated atom serials (heavy atoms)
    awk -v sel="$1" 'BEGIN { n = split(sel, s, ","); for (i = 1; i <= n; i++) { if (s[i] ~ /^[0-9]+-[0-9]+$/) { split(s[i], r, "-"); for (k = r[1]; k <= r[2]; k++) R[k] = 1 }
                                 else if (s[i] ~ /^[0-9]+$/) R[s[i] + 0] = 1; else N[s[i]] = 1 } }
        /^(ATOM|HETATM)/ { rn = substr($0, 18, 3); gsub(/ /, "", rn); ri = substr($0, 23, 4) + 0; el = substr($0, 77, 2); gsub(/ /, "", el)
                            an = substr($0, 13, 4); gsub(/ /, "", an); if (el == "") el = substr(an, 1, 1)
                            if ((ri in R || rn in N) && el != "H") { printf "%s%d", (c++ ? "," : ""), substr($0, 7, 5) + 0 } }' "$pdb"
}
L=$(atoms "$lig"); A=$(atoms "$sa"); [[ -n "$L" && -n "$A" ]] || die "selection matched no atoms"
echo "# PLUMED input written by plumed-input.sh ($cv), well-tempered metadynamics"
echo "lig: COM ATOMS=$L"
echo "pA:  COM ATOMS=$A"
if [[ $cv == funnel ]]; then
    [[ -n "$sb" ]] || die "funnel needs -b (axis point towards the solvent)"
    B=$(atoms "$sb"); [[ -n "$B" ]] || die "-b matched no atoms"
    cat <<EOT
pB:  COM ATOMS=$B
# fps.proj: position of the ligand along the axis, fps.ext: distance from the axis
fps: PROJECTION_ON_AXIS AXIS_ATOMS=pA,pB ATOM=lig
# allowed radius: cone (half-angle $alpha rad) below zcc = $zcc nm, cylinder of radius $rcyl nm above
rad: CUSTOM ARG=fps.proj VAR=z FUNC=step($zcc-z)*($rcyl+tan($alpha)*($zcc-z))+step(z-$zcc)*$rcyl PERIODIC=NO
exc: CUSTOM ARG=fps.ext,rad VAR=e,r FUNC=e-r PERIODIC=NO
funnel: UPPER_WALLS ARG=exc AT=0.0 KAPPA=20000.0
lwall:  LOWER_WALLS ARG=fps.proj AT=$zmin KAPPA=20000.0
uwall:  UPPER_WALLS ARG=fps.proj AT=$zmax KAPPA=20000.0
metad: METAD ARG=fps.proj SIGMA=$sigma HEIGHT=$height PACE=$pace BIASFACTOR=$gamma TEMP=$T GRID_MIN=$(calc "$zmin - 0.5") GRID_MAX=$(calc "$zmax + 0.5") FILE=HILLS
PRINT ARG=fps.proj,fps.ext,metad.bias,funnel.bias STRIDE=$pace FILE=COLVAR
EOT
else
    cat <<EOT
d: DISTANCE ATOMS=lig,pA
uwall: UPPER_WALLS ARG=d AT=$zmax KAPPA=20000.0
metad: METAD ARG=d SIGMA=$sigma HEIGHT=$height PACE=$pace BIASFACTOR=$gamma TEMP=$T GRID_MIN=0 GRID_MAX=$(calc "$zmax + 0.5") FILE=HILLS
PRINT ARG=d,metad.bias STRIDE=$pace FILE=COLVAR
EOT
fi
