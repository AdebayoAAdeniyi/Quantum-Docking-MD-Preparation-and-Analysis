#!/usr/bin/env bash
# =============================================================================
#  amber-analysis.sh - standard trajectory analysis with cpptraj, plus plots
#
#  usage: amber-analysis.sh -d SYSTEM_DIR [-y "prod.nc"] [-l :LIG] [-s 10] [-k 5]
#    -y  trajectories (default prod.nc; several: "equil.nc prod.nc")
#    -l  ligand mask (default :LIG or the first non-standard residue)
#    -s  analyse every s-th frame (default 1)
#    -k  number of clusters (k-means on ligand + pocket coordinates, default 5)
#  Produces in SYSTEM_DIR/analysis:
#    energy.dat      (from the mdout files: time, Etot, EPtot, Temp, Density)
#    rmsd-backbone.dat, rmsd-ligand.dat (ligand fitted on the protein backbone),
#    rmsf-ca.dat, radgyr.dat, hbonds-solute.dat (+ hbond-summary.dat),
#    ligand-contacts.dat (protein residues within 4 A of the ligand, % of frames),
#    clusters/ (representative structures), dssp-summary.dat
#    and PNG plots of all time series.
#  The old ptraj-* scripts each did one of these with ptraj; cpptraj does them
#  in a single pass over the trajectory.
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"
d=""; trajs="prod.nc"; lig=""; stride=1; k=5
while getopts "d:y:l:s:k:h" o; do
    case $o in d) d=$OPTARG ;; y) trajs=$OPTARG ;; l) lig=$OPTARG ;; s) stride=$OPTARG ;; k) k=$OPTARG ;; *) sed -n '2,20p' "$0"; exit 0 ;; esac
done
[[ -f "$d/solv.prmtop" ]] || die "no solv.prmtop in $d"
require cpptraj
a="$d/analysis"; mkdir -p "$a/clusters"
if [[ -z "$lig" ]]; then
    lig=$(awk '/%FLAG RESIDUE_LABEL/{f=1; getline; next} /%FLAG/{f=0} f' "$d/solv.prmtop" | tr -s ' ' '\n' |
          grep -vxE 'ALA|ARG|ASN|ASP|ASH|CYS|CYX|CYM|GLN|GLU|GLH|GLY|HIS|HIE|HID|HIP|ILE|LEU|LYS|LYN|MET|PHE|PRO|SER|THR|TRP|TYR|VAL|ACE|NME|NHE|WAT|HOH|Na\+|Cl-|K\+|' | head -1)
    [[ -n "$lig" ]] && lig=":$lig"
fi
haslig=0; [[ -n "$lig" ]] && haslig=1
# a protein is present when the topology has backbone residues
hasprot=0; prmtop_res() { awk '/%FLAG RESIDUE_LABEL/{f=1; getline; next} /%FLAG/{f=0} f' "$d/solv.prmtop" | tr -s ' ' '\n'; }
prmtop_res | grep -qxE 'ALA|GLY|LEU|SER|VAL|THR|LYS|ARG|ASP|GLU|PHE|TYR|TRP|HIS|HIE|HID|HIP|ILE|PRO|ASN|GLN|MET|CYS|CYX' && hasprot=1
prot='@CA,C,N'
{
    echo "parm $(realpath "$d/solv.prmtop")"
    for t in $trajs; do echo "trajin $(realpath "$d/$t") 1 last $stride"; done
    echo "autoimage"
    # RMSD on the protein backbone, or on the solute when there is no protein
    if ((hasprot)); then echo "rms Backbone $prot&!:WAT first out $a/rmsd-backbone.dat mass"
    else echo "rms Solute !:WAT,Na+,Cl-,K+&!@H= first out $a/rmsd-backbone.dat mass"; fi
    ((haslig)) && ((hasprot)) && echo "rms Ligand $lig&!@H= first nofit out $a/rmsd-ligand.dat"
    ((hasprot)) && echo "atomicfluct RMSF @CA out $a/rmsf-ca.dat byres"
    echo "radgyr Rg !:WAT,Na+,Cl-,K+&!@H= out $a/radgyr.dat mass"
    echo "hbond HB !:WAT,Na+,Cl-,K+ out $a/hbonds-solute.dat avgout $a/hbond-summary.dat"
    ((hasprot)) && echo "secstruct SS @CA out $a/dssp.dat sumout $a/dssp-summary.dat"
    if ((haslig)); then
        echo "nativecontacts $lig&!@H= \"!$lig&!:WAT,Na+,Cl-,K+&!@H=\" name Contacts distance 4.0 byresidue out $a/ligand-contacts-time.dat resout $a/ligand-contacts.dat"
        echo "cluster C1 kmeans clusters $k rms $lig&!@H= nofit sieve 5 summary $a/clusters/summary.dat info $a/clusters/info.dat repout $a/clusters/rep repfmt pdb"
    else
        echo "cluster C1 kmeans clusters $k rms !:WAT,Na+,Cl-,K+&!@H= sieve 5 summary $a/clusters/summary.dat info $a/clusters/info.dat repout $a/clusters/rep repfmt pdb"
    fi
    echo "run"
} > "$a/cpptraj.in"
cpptraj -i "$a/cpptraj.in" > "$a/cpptraj.log" 2>&1 || { tail -5 "$a/cpptraj.log" >&2; die "cpptraj failed ($a/cpptraj.log)"; }
# energies from every mdout in protocol order
for s in heat density equil prod; do [[ -f "$d/$s.out" ]] && awk -v st="$s" '
    /^ NSTEP =/ { t = $6 } / Etot   =/ { et = $3; ep = $9 } /TEMP\(K\) =/ { temp = $9 }
    / Density    =/ { den = $3 }
    /A V E R A G E S/ { exit }
    /^ -{20,}/ { if (t != "" && et != "") print st, t, et, ep, temp, (den == "" ? "nan" : den); t = et = den = "" }' "$d/$s.out"; done > "$a/energy.dat"
info "analysis in $a (ligand mask: ${lig:-none})"
"$(dirname "$0")/plot-md.sh" "$a" || true
