#!/usr/bin/env bash
# =============================================================================
#  vina-prepare.sh - receptor and ligand PDBQT files and the search box for
#                    AutoDock Vina
#
#  usage: vina-prepare.sh -r receptor.pdb [-x ref_ligand | -c "x y z" | -a "A:315:CA,A:381:CA"]
#                         [-p 8] [-s "20 20 20"] [-H 7.4] [-o vina-setup] [-l list -d dir | ligands...]
#    -r   receptor PDB with hydrogens (prepared, e.g. by PDBFixer, reduce or pdb4amber)
#    box (one of):
#    -x   reference ligand (crystal pose): box centre = its centre, size = its extent + 2p
#    -c   explicit centre "x y z" (with -s for the size)
#    -a   centre between two (or more) atoms "chain:resid:atom,...", the old
#         center-from-2-points recipe for pockets defined by two residues
#    -p   padding around the reference ligand in A (default 8)
#    -H   pH for ligand protonation with Open Babel (default 7.4; "none" keeps the input)
#  Ligands: any format, 2D or 3D; SMILES and flat structures get 3D coordinates.
#  Tools: Meeko (mk_prepare_receptor.py, mk_prepare_ligand.py; Forli lab), Open Babel.
#  Writes receptor.pdbqt, ligands/NAME.pdbqt, box.txt (Vina config) and box.pdb.
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"
rec=""; ref=""; centre=""; atoms=""; pad=8; size=""; ph=7.4; out=vina-setup; list=""; indir=.
while getopts "r:x:c:a:p:s:H:o:l:d:h" o; do
    case $o in r) rec=$OPTARG ;; x) ref=$OPTARG ;; c) centre=$OPTARG ;; a) atoms=$OPTARG ;; p) pad=$OPTARG ;; s) size=$OPTARG ;;
               H) ph=$OPTARG ;; o) out=$OPTARG ;; l) list=$OPTARG ;; d) indir=$OPTARG ;; *) sed -n '2,20p' "$0"; exit 0 ;; esac
done
shift $((OPTIND-1))
[[ -f "$rec" ]] || die "receptor PDB needed (-r)"
require obabel mk_prepare_ligand.py mk_prepare_receptor.py
mkdir -p "$out/ligands"
# ---- box ----
if [[ -n "$ref" ]]; then
    read -r cx cy cz sx sy sz <<< "$(xyz_body "$ref" | awk -v p="$pad" '
        NR==1 { for (k=2;k<=4;k++) { mn[k]=$k; mx[k]=$k } }
        { for (k=2;k<=4;k++) { if ($k<mn[k]) mn[k]=$k; if ($k>mx[k]) mx[k]=$k } }
        END { printf "%.3f %.3f %.3f %.1f %.1f %.1f", (mn[2]+mx[2])/2, (mn[3]+mx[3])/2, (mn[4]+mx[4])/2, mx[2]-mn[2]+2*p, mx[3]-mn[3]+2*p, mx[4]-mn[4]+2*p }')"
elif [[ -n "$atoms" ]]; then
    read -r cx cy cz <<< "$(awk -v sel="$atoms" 'BEGIN { n = split(sel, s, ","); for (i=1;i<=n;i++) { split(s[i], t, ":"); W[t[1] ":" t[2]+0 ":" t[3]] = 1 } }
        /^(ATOM|HETATM)/ { k = substr($0,22,1) ":" substr($0,23,4)+0 ":" substr($0,13,4); gsub(/ /, "", k)
                           if (k in W) { x += substr($0,31,8); y += substr($0,39,8); z += substr($0,47,8); m++ } }
        END { if (!m) exit 1; printf "%.3f %.3f %.3f", x/m, y/m, z/m }' "$rec")" || die "atoms $atoms not found in $rec"
    read -r sx sy sz <<< "${size:-22 22 22}"
elif [[ -n "$centre" ]]; then
    read -r cx cy cz <<< "$centre"; read -r sx sy sz <<< "${size:-22 22 22}"
else die "define the box with -x, -a or -c"; fi
[[ -n "$size" && -n "$ref" ]] && read -r sx sy sz <<< "$size"
printf 'center_x = %s\ncenter_y = %s\ncenter_z = %s\nsize_x = %s\nsize_y = %s\nsize_z = %s\n' "$cx" "$cy" "$cz" "$sx" "$sy" "$sz" > "$out/box.txt"
# box corners as a PDB for checking in a viewer
awk -v cx="$cx" -v cy="$cy" -v cz="$cz" -v sx="$sx" -v sy="$sy" -v sz="$sz" 'BEGIN { n = 0
    for (i=-1;i<=1;i+=2) for (j=-1;j<=1;j+=2) for (k=-1;k<=1;k+=2) printf "HETATM%5d  C   BOX X   1    %8.3f%8.3f%8.3f  1.00  0.00           C\n", ++n, cx+i*sx/2, cy+j*sy/2, cz+k*sz/2 }' > "$out/box.pdb"
info "box: centre ($cx, $cy, $cz), size $sx x $sy x $sz A"
# ---- receptor ----
mk_prepare_receptor.py -i "$rec" -o "$out/receptor" -p > "$out/receptor.log" 2>&1 ||
    mk_prepare_receptor.py --read_pdb "$rec" -o "$out/receptor" -p > "$out/receptor.log" 2>&1 || die "receptor preparation failed ($out/receptor.log)"
[[ -s "$out/receptor.pdbqt" ]] || die "no receptor.pdbqt (see $out/receptor.log)"
info "receptor: $out/receptor.pdbqt ($(grep -c '^ATOM\|^HETATM' "$out/receptor.pdbqt") atoms)"
# ---- ligands ----
mapfile -t ligs < <(gather_structures "$list" "$indir" "$@")
for f in "${ligs[@]}"; do
    n=$(strip_ext "$f"); tmp="$out/ligands/$n.sdf"
    if [[ $ph == none ]]; then to_format "$f" "$tmp"; else to_format "$f" "$tmp" -p "$ph"; fi
    if mk_prepare_ligand.py -i "$tmp" -o "$out/ligands/$n.pdbqt" > "$out/ligands/$n.log" 2>&1; then
        info "ligand $n: $(grep -c '^ATOM\|^HETATM' "$out/ligands/$n.pdbqt") atoms, $(awk '/^TORSDOF/{print $2}' "$out/ligands/$n.pdbqt") rotatable bonds"
    else warn "ligand $n failed (see $out/ligands/$n.log)"; fi
done
