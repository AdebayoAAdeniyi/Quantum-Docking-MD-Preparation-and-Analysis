#!/usr/bin/env bash
# =============================================================================
#  multiwfn-batch.sh - run Multiwfn non-interactively over a set of molecules
#
#  Multiwfn (Tian Lu, http://sobereva.com/multiwfn; T. Lu and F. Chen,
#  J. Comput. Chem. 2012, 33, 580) is menu driven. Each task below is the menu
#  key sequence that the original workflow used, fed on standard input; results
#  are renamed after the molecule and the numbers collected in one CSV.
#  The menu numbers follow Multiwfn 3.8. If a newer version moves a menu,
#  record the new key sequence once with -k.
#
#  usage: multiwfn-batch.sh -t TASK [-l list] [-p parent] [-o outdir] [-q 2|3]
#   TASK         input (in parent/NAME/)                  output
#   fukui        NAME-fukui-{N,Nplus1,Nminus1}.fchk       f+ and f- density cubes
#                (made by gaussian-input.sh -F)
#   mesp         NAME.fchk                                ESP on the vdW surface: Vmin, Vmax,
#                                                         vertices and extrema as pdb
#   homo-lumo    NAME.fchk                                HOMO and LUMO cubes
#   overlap      NAME.fchk                                HOMO/LUMO centroid distance and overlap
#   hole-electron NAME-td.fchk + NAME.log (TD job)        hole/electron cubes, centroid distance,
#                                                         exciton Coulomb energy  [-s state]
#   custom       NAME.fchk + key file -k FILE             anything; "NAME" in the key file is
#                                                         replaced by the molecule name
#   -q  grid quality for cubes: 2 medium, 3 high (default 3)
#  The isosurfaces can be drawn with VMD; Tian Lu's orbitalshow.vmd script from
#  the Multiwfn manual is the convenient way and is not redistributed here.
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"
task=""; list=""; parent=gaussian-jobs; out=multiwfn-results; keys=""; state=1; grid=3
while getopts "t:l:p:o:k:s:q:h" o; do
    case $o in t) task=$OPTARG ;; l) list=$OPTARG ;; p) parent=$OPTARG ;; o) out=$OPTARG ;; k) keys=$OPTARG ;;
               s) state=$OPTARG ;; q) grid=$OPTARG ;; *) sed -n '2,29p' "$0"; exit 0 ;; esac
done
[[ -n "$task" ]] || die "choose -t (see -h)"
MW=${MULTIWFN_EXE:-Multiwfn}; require "$MW"
if [[ -n "$list" ]]; then mapfile -t names < <(read_list "$list")
else mapfile -t names < <(find "$parent" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort); fi
mkdir -p "$out"; out=$(cd "$out" && pwd)
csv="$out/$task.csv"
case $task in
    mesp)          echo "name,Vmin_kcal_mol,Vmax_kcal_mol" > "$csv" ;;
    overlap)       echo "name,HOMO,LUMO,centroid_distance_A,overlap_norm,overlap_square" > "$csv" ;;
    hole-electron) echo "name,state,E_exc_eV,D_index_A,Coulomb_attraction_eV" > "$csv" ;;
    *) : > "$csv" ;;
esac
homo_index() {   # number of the alpha HOMO, read from Multiwfn's first screen
    printf '%s\n0\n' "$1" | "$MW" 2>/dev/null |
        awk '/is HOMO of alpha spin/ {print $3; exit} /HOMO, energy/ && /Orbital/ {print $3; exit}' | tr -dc '0-9'
}
for n in "${names[@]}"; do
    d="$parent/$n"; w="$out/$n"; mkdir -p "$w"
    case $task in
    fukui)
        for s in N Nplus1 Nminus1; do [[ -f "$d/$n-fukui-$s.fchk" ]] || { warn "$n: $n-fukui-$s.fchk missing"; continue 2; }; done
        # f+ = rho(N+1) - rho(N) ; f- = rho(N) - rho(N-1)   (custom operation "-" on the electron density)
        (cd "$w" && printf '%s\n5\n0\n1\n-,%s\n1\n%s\n2\n' "$d/$n-fukui-Nplus1.fchk" "$d/$n-fukui-N.fchk" "$grid" | "$MW" > fplus.out 2>&1 && mv density.cub "$n-fplus.cub")
        (cd "$w" && printf '%s\n5\n0\n1\n-,%s\n1\n%s\n2\n' "$d/$n-fukui-N.fchk" "$d/$n-fukui-Nminus1.fchk" "$grid" | "$MW" > fminus.out 2>&1 && mv density.cub "$n-fminus.cub")
        info "$n: $w/$n-fplus.cub, $n-fminus.cub" ;;
    mesp)
        [[ -f "$d/$n.fchk" ]] || { warn "$n: $n.fchk missing"; continue; }
        # 12 surface analysis, 0 start (ESP), 6 vertices pdb, 2 extrema pdb, 1 extrema txt, 5 molecule pdb
        (cd "$w" && printf '%s\n12\n0\n6\n2\n1\n5\n%s-geom.pdb\n' "$d/$n.fchk" "$n" | "$MW" > mesp.out 2>&1)
        [[ -f "$w/vtx.pdb" ]] && mv "$w/vtx.pdb" "$w/$n-mesp-vertices.pdb"
        [[ -f "$w/surfanalysis.pdb" ]] && mv "$w/surfanalysis.pdb" "$w/$n-extrema.pdb"
        # global minimum and maximum of the ESP on the surface (kcal/mol column of surfanalysis.txt)
        if [[ -f "$w/surfanalysis.txt" ]]; then
            awk -v n="$n" '$1 ~ /^\*?[0-9]+$/ && NF >= 5 { v = $4 + 0; if (min == "" || v < min) min = v; if (max == "" || v > max) max = v }
                 END { print n "," min "," max }' "$w/surfanalysis.txt" >> "$csv"
        fi
        info "$n: ESP surface analysis in $w" ;;
    homo-lumo)
        h=$(homo_index "$d/$n.fchk"); [[ -n "$h" ]] || { warn "$n: HOMO not found"; continue; }
        (cd "$w" && printf '%s\n200\n3\n%s-%s\n%s\n1\n' "$d/$n.fchk" "$h" "$((h + 1))" "$grid" | "$MW" > homo-lumo.out 2>&1)
        info "$n: orbitals $h (HOMO) and $((h + 1)) (LUMO) written as orb*.cub in $w" ;;
    overlap)
        h=$(homo_index "$d/$n.fchk"); [[ -n "$h" ]] || { warn "$n: HOMO not found"; continue; }
        (cd "$w" && printf '%s\n100\n11\n%s,%s\n0,0\n' "$d/$n.fchk" "$h" "$((h + 1))" | "$MW" > overlap.out 2>&1)
        awk -v n="$n" -v h="$h" '/Centroid distance between the two orbitals:/ {d = $(NF-1)}
             /Overlap integral of norm of the two orbitals:/ {on = $NF} /Overlap integral of square of the two orbitals:/ {os = $NF}
             END { print n "," h "," h + 1 "," d "," on "," os }' "$w/overlap.out" >> "$csv" ;;
    hole-electron)
        f="$d/$n-td.fchk"; [[ -f "$f" ]] || f="$d/$n.fchk"; lg="$d/$n.log"
        [[ -f "$f" && -f "$lg" ]] || { warn "$n: need the TD fchk and log"; continue; }
        # 18 excitation analysis, 1 hole-electron, log file, state, 1 visualise, grid, 10/11 export hole/electron, 13 transition density
        (cd "$w" && printf '%s\n18\n1\n%s\n%s\n1\n%s\n10\n1\n11\n1\n13\n0\n0\nq\n' "$f" "$lg" "$state" "$grid" | "$MW" > hole-electron.out 2>&1)
        for c in hole electron transdens; do [[ -f "$w/$c.cub" ]] && mv "$w/$c.cub" "$w/$n-S$state-$c.cub"; done
        awk -v n="$n" -v s="$state" '/Excitation energy of this state:/ {e = $(NF-1)} /D index:/ {di = $(NF-1)}
             /Coulomb attractive energy:/ {c = $(NF-2)} END { print n "," s "," e "," di "," c }' "$w/hole-electron.out" >> "$csv" ;;
    custom)
        [[ -f "$keys" ]] || die "custom needs a key file (-k)"
        (cd "$w" && sed "s|NAME|$n|g; s|PARENT|$d|g" "$keys" | "$MW" > custom.out 2>&1)
        info "$n: $w/custom.out" ;;
    *) die "unknown task $task" ;;
    esac
done
[[ -s "$csv" ]] && info "table: $csv" || rm -f "$csv"
