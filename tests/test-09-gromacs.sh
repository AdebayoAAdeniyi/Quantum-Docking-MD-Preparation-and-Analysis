#!/usr/bin/env bash
# GROMACS: convert an Amber ligand box, build the index groups and mdp files,
# and run a real energy minimisation.
source "$(dirname "$0")/testlib.sh"
need gmx; need tleap; need antechamber
G=$QDM_ROOT/09-gromacs; A=$QDM_ROOT/08-amber
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT; cd "$W"
fetch "https://raw.githubusercontent.com/ccsb-scripps/AutoDock-Vina/develop/example/basic_docking/data/1iep_ligand.sdf" 1iep_ligand.sdf || exit 0
cp "$DATA/1iep_ligand.sdf" STI.sdf
"$A/prepare/amber-ligand-params.sh" -m gas -o lig STI.sdf 2>/dev/null
"$A/prepare/amber-build.sh" -n lw -L lig/STI -b 8 -o asys 2>/dev/null
check "Amber system built"   "$([[ -s asys/solv.prmtop ]] && echo yes)" yes
"$G/prepare/gmx-build.sh" -n lw -a asys -o gsys 2>/dev/null
check "GROMACS topology"     "$([[ -s gsys/topol.top ]] && echo yes)" yes
check "same atom count"      "$(sed -n 2p gsys/system.gro | tr -d ' ')" "$(awk '/%FLAG POINTERS/{getline; getline; print $1; exit}' asys/solv.prmtop)"
check "LIG group"            "$(grep -c "^\[ LIG \]" gsys/index.ndx)" 1
check "Water_and_ions group" "$(grep -c '^\[ Water_and_ions \]' gsys/index.ndx)" 1
"$G/run/gmx-md-inputs.sh" -d gsys -T 300 -V 2 -N 2 -P 0.002 -f 1 -A 2>/dev/null
check "four mdp files"       "$(ls gsys/{em,nvt,npt,md}.mdp 2>/dev/null | wc -l)" 4
check "C-rescale barostat"   "$(grep -c 'pcoupl = C-rescale' gsys/npt.mdp)" 1
check "annealing in nvt"     "$(grep -c '^annealing' gsys/nvt.mdp)" 4
check "no annealing in md"   "$(grep -c '^annealing' gsys/md.mdp)" 0
check "position restraints"  "$(ls gsys/posre-*.itp | wc -l | awk '{print ($1>0)?"yes":"no"}')" yes
check "POSRES include"       "$(grep -c 'ifdef POSRES' gsys/topol.top | awk '{print ($1>0)?"yes":"no"}')" yes
(cd gsys && gmx grompp -f em.mdp -c system.gro -r system.gro -p topol.top -n index.ndx -o em.tpr -maxwarn 1 > grompp.log 2>&1 &&
          gmx mdrun -deffnm em -nt 1 -nsteps 50 > mdrun.log 2>&1)
check "minimisation runs"    "$([[ -s gsys/em.gro ]] && echo yes)" yes
finish
