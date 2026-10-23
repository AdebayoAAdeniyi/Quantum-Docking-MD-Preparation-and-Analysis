#!/usr/bin/env bash
# PLUMED: a real well-tempered metadynamics run on a double well whose barrier is
# known analytically (5 energy units), plus the funnel input generator.
source "$(dirname "$0")/testlib.sh"
P=$QDM_ROOT/12-plumed-metadynamics
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT; cd "$W"
echo "== funnel input generation"
cat > c.pdb <<'EOT'
ATOM      1  CA  ALA A 315      10.000  10.000  10.000  1.00  0.00           C
ATOM      2  CB  ALA A 315      11.000  10.000  10.000  1.00  0.00           C
ATOM      3  CA  ALA A  40      10.000  10.000  25.000  1.00  0.00           C
ATOM      4  C1  LIG A 999      10.500  10.500  12.000  1.00  0.00           C
ATOM      5  N1  LIG A 999      11.500  10.500  12.000  1.00  0.00           N
ATOM      6  H1  LIG A 999      12.000  10.500  12.000  1.00  0.00           H
EOT
"$P/plumed-input.sh" -p c.pdb -l LIG -a 315 -b 40 -z 1.8 -r 0.1 > plumed-funnel.dat 2>/dev/null
check "ligand COM, hydrogens excluded" "$(grep -c 'lig: COM ATOMS=4,5' plumed-funnel.dat)" 1
check "projection on the axis"   "$(grep -c 'PROJECTION_ON_AXIS' plumed-funnel.dat)" 1
check "funnel wall"              "$(grep -c 'funnel: UPPER_WALLS' plumed-funnel.dat)" 1
check "well-tempered metad"      "$(grep -c 'BIASFACTOR=10' plumed-funnel.dat)" 1
"$P/plumed-input.sh" -p c.pdb -l LIG -a 315 -c distance -M 3.0 > plumed-dist.dat 2>/dev/null
check "distance CV alternative"  "$(grep -c '^d: DISTANCE' plumed-dist.dat)" 1
need plumed
echo "== metadynamics on a double well (analytic barrier = 5)"
cp "$QDM_ROOT/examples/metad-double-well/in.pesmd" "$QDM_ROOT/examples/metad-double-well/plumed.dat" .
plumed pesmd < in.pesmd > pesmd.log 2>&1
check "hills deposited"          "$(grep -vc '^#' HILLS | awk '{print ($1>500)?"yes":"no"}')" yes
"$P/plumed-fes.sh" -i HILLS -k 1.0 -n 4 -o fes 2>/dev/null
b=$(awk '!/^#/ && $1 > -0.05 && $1 < 0.05 {print $2; exit}' fes.dat)
m=$(awk '!/^#/ && $1 > 0.95 && $1 < 1.05 {print $2; exit}' fes.dat)
check "barrier at x=0 is 5"      "$b" 5.0 0.8
check "minimum at x=1 is 0"      "$m" 0.0 0.3
check "convergence series"       "$(grep -c '^$' fes-convergence.dat | awk '{print ($1>3)?"yes":"no"}')" yes
finish
