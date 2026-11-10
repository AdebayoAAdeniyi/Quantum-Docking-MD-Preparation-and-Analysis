#!/usr/bin/env bash
# Structure preparation: everything here runs, because Open Babel is free.
#   build-3d-from-smiles.sh  3D generation and conformer search from SMILES
#   convert-structures.sh    format round trip, hydrogens, pH, centring
#   ligand-pdb-for-md.sh     unique atom names and one residue, as tleap needs
#   split-multistructure.sh  one file per molecule, nothing lost or overwritten
source "$(dirname "$0")/testlib.sh"
need obabel
S=$QDM_ROOT/00-structure-preparation
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT; cd "$W"

echo "== build-3d-from-smiles.sh"
printf '1-H\tOc1ccccc1C(=O)c1ccccc1\nbenzene\tc1ccccc1\n' > smiles.txt
"$S/build-3d-from-smiles.sh" -i smiles.txt -o built -f xyz,pdb -n 5 -s 200 >/dev/null 2>&1
check "both molecules built"      "$(ls built/*.xyz 2>/dev/null | wc -l)" 2
check "benzene has 12 atoms"      "$(head -1 built/benzene.xyz 2>/dev/null | tr -d ' \r')" 12
# 2-hydroxybenzophenone is C13H10O2: 25 atoms
check "1-H has 25 atoms"          "$(head -1 built/1-H.xyz 2>/dev/null | tr -d ' \r')" 25
check "coordinates are 3D"        "$(awk 'NR==3{print ($3!=0)||($4!=0)?1:0}' built/benzene.xyz)" 1

echo "== convert-structures.sh"
# a round trip must not lose or invent atoms
"$S/convert-structures.sh" -o conv -f pdb,mol2 built/benzene.xyz >/dev/null 2>&1
check "pdb written"               "$(ls conv/benzene.pdb 2>/dev/null | wc -l)" 1
check "mol2 written"              "$(ls conv/benzene.mol2 2>/dev/null | wc -l)" 1
check "atom count survives"       "$(grep -c '^\(ATOM\|HETATM\)' conv/benzene.pdb)" 12
check "carbons survive"           "$(obabel conv/benzene.pdb -oxyz 2>/dev/null | awk '$1=="C"{n++}END{print n}')" 6
# centring must put the centroid at the origin
"$S/convert-structures.sh" -o cent -f xyz -c built/benzene.xyz >/dev/null 2>&1
check "centred at the origin"     "$(awk 'NR>2{x+=$2;y+=$3;z+=$4;n++}END{printf "%.3f", (x*x+y*y+z*z)/n}' cent/benzene.xyz)" 0.000 0.01
# pH protonation must change the hydrogen count of an acid
obabel -:"CC(=O)O" -oxyz -O acid.xyz --gen3d >/dev/null 2>&1
"$S/convert-structures.sh" -o ph2  -f xyz -p 2  acid.xyz >/dev/null 2>&1
"$S/convert-structures.sh" -o ph11 -f xyz -p 11 acid.xyz >/dev/null 2>&1
h2=$(awk 'NR>2 && $1=="H"{n++}END{print n+0}' ph2/acid.xyz)
h11=$(awk 'NR>2 && $1=="H"{n++}END{print n+0}' ph11/acid.xyz)
check "acetic acid loses its proton above the pKa" "$((h2 > h11 ? 1 : 0))" 1

echo "== ligand-pdb-for-md.sh"
"$S/ligand-pdb-for-md.sh" built/1-H.xyz LIG > lig.pdb 2>/dev/null
n_at=$(grep -c '^\(ATOM\|HETATM\)' lig.pdb)
n_uniq=$(awk '/^(ATOM|HETATM)/{print substr($0,13,4)}' lig.pdb | sort -u | wc -l)
check "every atom kept"           "$n_at" 25
check "every atom name is unique" "$n_uniq" "$n_at"
check "one residue name"          "$(awk '/^(ATOM|HETATM)/{r=substr($0,18,3); gsub(/ /,"",r); print r}' lig.pdb | sort -u | tr -d '\n')" LIG
check "one residue number"        "$(awk '/^(ATOM|HETATM)/{print substr($0,23,4)+0}' lig.pdb | sort -u | wc -l)" 1
check "names start with the element" \
    "$(awk '/^(ATOM|HETATM)/{an=substr($0,13,4); gsub(/ /,"",an); el=substr($0,77,2); gsub(/ /,"",el);
            if (el != "" && substr(an,1,length(el)) != el) bad++} END{print bad+0}' lig.pdb)" 0

echo "== split-multistructure.sh"
# three named structures in one sdf, with one title repeated
for n in poseA poseB poseA; do obabel -:"c1ccccc1" -osdf --gen3d --title "$n" >> multi.sdf 2>/dev/null; done
"$S/split-multistructure.sh" multi.sdf -o split >/dev/null 2>&1
check "nothing lost in the split" "$(ls split/*.sdf 2>/dev/null | wc -l)" 3
check "repeated title not overwritten" "$(ls split/ | grep -c 'poseA')" 2
"$S/split-multistructure.sh" multi.sdf -o best -k 1 >/dev/null 2>&1
check "-k 1 keeps one per title"  "$(ls best/*.sdf 2>/dev/null | wc -l)" 2

finish
