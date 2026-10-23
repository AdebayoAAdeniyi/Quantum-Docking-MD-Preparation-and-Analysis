#!/usr/bin/env bash
# =============================================================================
#  gmx-md-inputs.sh - mdp files and a restartable job script for the usual
#                     GROMACS protocol
#
#   em     steepest descent minimisation
#   nvt    heating at constant volume, solute heavy atoms restrained; with -A a
#          simulated annealing ramp 0 -> T K (the original workflow heated this way)
#   npt    density equilibration, restraints kept, C-rescale barostat
#   md     production, no restraints, V-rescale + C-rescale, 2 fs, LINCS on H bonds
#
#  usage: gmx-md-inputs.sh -d SYSTEM_DIR [-T 300] [-V 100] [-N 100] [-P 10] [-f 10] [-A] [-g]
#    -V/-N   NVT / NPT length (ps)    -P  production (ns)   -f  frame interval (ps)
#    -A      annealing during NVT     -g  offload to the GPU (-nb gpu -pme gpu -update gpu)
#  Temperature coupling groups: Protein_LIG (or Protein) and Water_and_ions.
#  Restraints: heavy atoms of the solute (posre.itp made with gmx genrestr and
#  included under #ifdef POSRES when the topology has none).
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"
source "$QDM_ROOT/lib/jobscript.sh"; load_cluster_conf
d=""; T=300; nvt=100; npt=100; prod=10; frame=10; anneal=0; gpu=0
while getopts "d:T:V:N:P:f:Agh" o; do
    case $o in d) d=$OPTARG ;; T) T=$OPTARG ;; V) nvt=$OPTARG ;; N) npt=$OPTARG ;; P) prod=$OPTARG ;; f) frame=$OPTARG ;;
               A) anneal=1 ;; g) gpu=1 ;; *) sed -n '2,19p' "$0"; exit 0 ;; esac
done
[[ -f "$d/topol.top" && -f "$d/system.gro" ]] || die "no topol.top/system.gro in $d"
GMX=${GMX:-gmx}
grp=Protein; grep -q "^\[ Protein_LIG \]" "$d/index.ndx" 2>/dev/null && grp=Protein_LIG
grep -q "^\[ Protein \]" "$d/index.ndx" 2>/dev/null || grp=non-Water
sol=Water_and_ions; grep -q "^\[ Water_and_ions \]" "$d/index.ndx" 2>/dev/null || sol=Water
st() { awk -v ps="$1" 'BEGIN{printf "%d", ps / 0.002}'; }
nst=$(st "$frame")
common="cutoff-scheme = Verlet
nstlist = 20
coulombtype = PME
rcoulomb = 1.0
rvdw = 1.0
DispCorr = EnerPres
constraints = h-bonds
constraint-algorithm = lincs"
printf 'integrator = steep\nemtol = 1000.0\nemstep = 0.01\nnsteps = 50000\n%s\n' "$common" > "$d/em.mdp"
{
    echo "define = -DPOSRES"
    echo "integrator = md"; echo "dt = 0.002"; echo "nsteps = $(st "$nvt")"
    echo "nstxout-compressed = $nst"; echo "nstenergy = 500"; echo "nstlog = 500"
    echo "$common"
    echo "tcoupl = V-rescale"; echo "tc-grps = $grp $sol"; echo "tau-t = 0.1 0.1"; echo "ref-t = $T $T"
    echo "pcoupl = no"
    echo "gen-vel = yes"; echo "gen-temp = $( ((anneal)) && echo 0 || echo "$T" )"; echo "gen-seed = -1"
    echo "continuation = no"
    if ((anneal)); then
        a=$(calc "$nvt * 0.8")
        echo "annealing = single single"; echo "annealing-npoints = 2 2"
        echo "annealing-time = 0 $a 0 $a"; echo "annealing-temp = 0 $T 0 $T"
    fi
} > "$d/nvt.mdp"
sed -e 's/^pcoupl = no/pcoupl = C-rescale\npcoupltype = isotropic\ntau-p = 2.0\nref-p = 1.0\ncompressibility = 4.5e-5\nrefcoord-scaling = com/' \
    -e 's/^gen-vel = yes/gen-vel = no/' -e '/^gen-temp/d' -e '/^gen-seed/d' -e 's/^continuation = no/continuation = yes/' \
    -e "s/^nsteps = .*/nsteps = $(st "$npt")/" -e '/^annealing/d' "$d/nvt.mdp" > "$d/npt.mdp"
sed -e '/^define/d' -e "s/^nsteps = .*/nsteps = $(st "$(calc "$prod * 1000")")/" -e '/^refcoord-scaling/d' "$d/npt.mdp" > "$d/md.mdp"

# position restraints for topologies without them (e.g. converted from Amber):
# heavy atoms of every molecule type that is not water or an ion
if ! grep -q "POSRES" "$d/topol.top"; then
    awk -v dir="$d" '
        function close_mol() { if (inm && nr) { print "#ifdef POSRES"; print "#include \"posre-" cur ".itp\""; print "#endif"; print "" } inm = 0; nr = 0 }
        /^\[ *moleculetype *\]/ { close_mol(); print; getline; print; while ($0 ~ /^;/) { getline; print }
                                   cur = $1; inm = (cur !~ /^(WAT|SOL|HOH|Na\+|Cl-|NA|CL|K\+|Na|Cl|K)$/); next }
        /^\[ *atoms *\]/ { inat = inm; print; next }
        /^\[/ { inat = 0 }
        /^\[ *system *\]/ { close_mol() }
        inat && $1 ~ /^[0-9]+$/ && $5 !~ /^H/ { f = dir "/posre-" cur ".itp"
            if (!(f in seen)) { seen[f] = 1; print "[ position_restraints ]" > f }
            printf "%6d 1 1000 1000 1000\n", $1 > f; nr++ }
        { print }' "$d/topol.top" > "$d/topol.tmp" && mv "$d/topol.tmp" "$d/topol.top"
    n=$(ls "$d"/posre-*.itp 2>/dev/null | wc -l)
    ((n)) && info "position restraints written for $n solute molecule type(s)" || warn "no solute molecule found for restraints"
fi
gpuopt=""; ((gpu)) && gpuopt="-nb gpu -pme gpu -update gpu"
cmds=("GMX=\"$GMX\"" "NDX=\$( [[ -f index.ndx ]] && echo '-n index.ndx' )"
      "step() { local s=\$1 c=\$2; [[ -s \$s.gro ]] && { echo \"\$s done\"; return 0; }"
      "  \$GMX grompp -f \$s.mdp -c \$c -r \$c -p topol.top \$NDX -o \$s.tpr -maxwarn 1 && \$GMX mdrun -deffnm \$s -nt \$NCORES $gpuopt -cpi \$s.cpt; }"
      "step em system.gro && step nvt em.gro && step npt nvt.gro && step md npt.gro")
write_jobscript "gmx-$(basename "$d")" "$PWD/$d" GROMACS_SETUP "${cmds[@]}" > "$d/run-md.sh"; chmod +x "$d/run-md.sh"
info "$d: em, nvt ($nvt ps$( ((anneal)) && echo ', annealing')), npt ($npt ps), md ($prod ns); coupling groups $grp / $sol"
