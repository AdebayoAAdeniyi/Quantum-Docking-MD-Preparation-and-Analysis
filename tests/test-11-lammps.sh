#!/usr/bin/env bash
# LAMMPS: a short real run of the bead-spring melt, and the analysis scripts
# checked against the stored results of the full run in examples/.
source "$(dirname "$0")/testlib.sh"
L=$QDM_ROOT/11-lammps; E=$QDM_ROOT/examples/lammps-bead-spring
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT; cd "$W"
echo "== analysis on the stored example run"
"$L/analysis/lammps-tg.sh" -i "$E/tg-density.dat" -o tg 2>/dev/null
check "Tg from 9 temperatures"  "$(awk -F, 'NR==2{print $1}' tg.csv)" 0.5692 0.001
check "melt expands more than glass" "$(awk -F, 'NR==2{print ($7 > $6) ? "yes" : "no"}' tg.csv)" yes
"$L/analysis/lammps-viscosity.sh" -i "$E/viscosity-running.dat" -o eta 2>/dev/null
check "Green-Kubo plateau"      "$(awk -F, 'NR==2{print $1}' eta.csv)" 3.2263 0.01
"$L/analysis/lammps-elastic.sh" -i "$E/stress-strain-T0.3.dat" -m 0.03 -o el 2>/dev/null
check "Young's modulus"         "$(awk -F, 'NR==2{print $1}' el.csv)" 7.05 0.1
echo "== build and a short real run"
"$L/prepare/build-bead-spring.sh" -c 10 -n 10 -r 0.85 -o poly 2>/dev/null
check "100 beads"               "$(awk '/atoms/{print $1}' poly/system.data)" 100
check "90 bonds"                "$(awk '/bonds/{print $1}' poly/system.data)" 90
"$L/prepare/lammps-inputs.sh" -s poly/system.in -D poly/system.data -u lj -x 0.02 -o run 2>/dev/null
check "five inputs + system.in"  "$(ls run/*.in | grep -cv params)" 6
check "params.in has T"         "$(awk '$2=="T"{print $4}' run/params.in)" 1.0
need lmp
(cd run && lmp -in equilibrate.in -log eq.log > eq.out 2>&1)
check "equilibration finished"  "$([[ -s run/equilibrated.data ]] && echo yes)" yes
check "density is physical"     "$(tail -1 run/equilibrate-density.dat | awk '{print $2}')" 0.85 0.25
finish
