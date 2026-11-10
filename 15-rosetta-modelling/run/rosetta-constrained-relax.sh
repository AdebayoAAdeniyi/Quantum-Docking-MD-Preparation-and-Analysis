#!/usr/bin/env bash
# =============================================================================
#  rosetta-constrained-relax.sh - relax under experimental distance restraints
#
#  FastRelax with AtomPair restraints, for imposing contacts you have evidence
#  for: a crosslink, an NOE, a salt bridge or hydrogen bond from mutagenesis, a
#  spin label distance. This is the step that turns a plausible assembly into
#  one consistent with the measurements.
#
#  usage: rosetta-constrained-relax.sh -p COMPLEX.pdb -c RESTRAINTS [options]
#    -p FILE    the structure to relax
#    -c FILE    restraint file, Rosetta constraint syntax or the short form below
#    -o DIR     output directory (default: constrained_relax)
#    -n N       number of structures (default: 100)
#    -W W       constraint weight (default: 4)
#    -f         FastRelax instead of dual-space relax (faster, less thorough)
#    -q         write the job script and stop, for the queue
#    -h         this help
#
#  The restraint file
#  ------------------
#  Rosetta's own AtomPair syntax is accepted as given:
#
#    AtomPair  H  163A  H  167B  BOUNDED  1.5  6.5  1.0  NOE
#
#  A short form is also accepted, one restraint per line, and is expanded:
#
#    163A:H  167B:H  1.5  6.5      # lower and upper bound in Angstrom
#    176A:ND2  537B:CG  2.6  3.4
#
#  Residue numbering with a chain letter, as in 163A, is PDB numbering. Rosetta
#  accepts that form, and it is what you want for a multi-chain assembly, since
#  pose numbering runs continuously across chains and the residue you mean is
#  not the residue you typed. The script checks every atom named in the file
#  actually exists in the structure, which is the failure the original script
#  could not catch: a restraint naming an atom that is not there is dropped
#  silently, and the run completes having enforced nothing.
#
#  Weight
#  ------
#  The default of 4 follows the Rosetta documentation's suggestion for folding
#  under restraints. It is strong. A weight that forces a contact the structure
#  cannot accommodate will distort the rest of the fold to get there, so check
#  the restraint energies in the score table afterwards and not only the total.
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../../lib/common.sh"
source "$QDM_ROOT/config/cluster.conf"
[[ -f "$QDM_ROOT/config/cluster.local.conf" ]] && source "$QDM_ROOT/config/cluster.local.conf"

pdb=""; cst=""; out="constrained_relax"; nstruct=100; weight=4; fast=0; queue=0
while getopts "p:c:o:n:W:fqh" o; do case $o in
    p) pdb=$OPTARG ;; c) cst=$OPTARG ;; o) out=$OPTARG ;; n) nstruct=$OPTARG ;;
    W) weight=$OPTARG ;; f) fast=1 ;; q) queue=1 ;;
    h) sed -n '2,/^# ====/p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) die "unknown option; -h for usage" ;;
esac; done

[[ -n "$pdb" ]] || die "a structure is required (-p); -h for usage"
[[ -f "$pdb" ]] || die "not found: $pdb"
[[ -n "$cst" ]] || die "a restraint file is required (-c); -h for usage"
[[ -f "$cst" ]] || die "restraint file not found: $cst"

mkdir -p "$out/input_files" "$out/output_files"
base=$(strip_ext "$pdb")
cp "$pdb" "$out/input_files/$base.pdb"

# ---------------------------------------------------------------------------
# Expand the short form into Rosetta AtomPair lines
# ---------------------------------------------------------------------------
awk '
    /^[[:space:]]*#/ || /^[[:space:]]*$/ { print; next }
    $1 == "AtomPair" || $1 == "Dihedral" || $1 == "Angle" || $1 == "CoordinateConstraint" { print; next }
    NF >= 4 && $1 ~ /:/ && $2 ~ /:/ {
        split($1, a, ":"); split($2, b, ":")
        printf "AtomPair %6s %8s %6s %8s  BOUNDED %6s %6s %4s  experimental\n",
               a[2], a[1], b[2], b[1], $3, $4, (NF >= 5 ? $5 : "1.0")
        next
    }
    { printf "line %d of the restraint file was not understood: %s\n", NR, $0 > "/dev/stderr"; bad = 1 }
    END { exit bad }
' "$cst" > "$out/input_files/$base.cst" || die "fix the restraint file and run again"

n_cst=$(grep -c '^AtomPair' "$out/input_files/$base.cst" || true)
(( n_cst > 0 )) || die "no restraints found in $cst"
info "$n_cst restraint(s) prepared"

# ---------------------------------------------------------------------------
# Check every restrained atom exists. Rosetta drops a restraint that names a
# missing atom without stopping, so the run "succeeds" having enforced nothing.
# ---------------------------------------------------------------------------
awk -v pdbfile="$out/input_files/$base.pdb" '
    BEGIN {
        while ((getline line < pdbfile) > 0) {
            if (line !~ /^(ATOM|HETATM)/) continue
            an = substr(line, 13, 4); gsub(/ /, "", an)
            rn = substr(line, 23, 4) + 0
            ch = substr(line, 22, 1)
            have[rn ch ":" an] = 1
            have[rn ":" an] = 1
        }
        close(pdbfile)
    }
    $1 == "AtomPair" {
        split($3, r1, ""); k1 = $3 ":" $2
        k2 = $5 ":" $4
        if (!(k1 in have)) { printf "  restraint %d: atom %s of residue %s not found\n", ++n, $2, $3 > "/dev/stderr"; bad = 1 }
        if (!(k2 in have)) { printf "  restraint %d: atom %s of residue %s not found\n", ++n, $4, $5 > "/dev/stderr"; bad = 1 }
    }
    END { exit bad }
' "$out/input_files/$base.cst" || die "the restraints above name atoms that are not in $pdb
Backbone amide hydrogen is H in Rosetta full atom naming. Check the chain
letters and that hydrogens are present: a structure straight from a crystal
has none, so a hydrogen restraint cannot be satisfied."
info "every restrained atom is present in the structure"

{
    cat <<EOF
# constrained relax, generated by rosetta-constrained-relax.sh
-in:file:s input_files/$base.pdb
-nstruct $nstruct

EOF
    if ((fast)); then
        printf -- '-relax:fast\n'
    else
        printf -- '-relax:dualspace\n-relax:minimize_bond_angles\n'
    fi
    cat <<EOF
-set_weights cart_bonded 0.5 pro_close 0
-default_max_cycles 200

# restraints. Centroid and full atom stages take separate files; both are
# given here so the restraints apply throughout rather than only at the end.
-constraints:cst_file input_files/$base.cst
-constraints:cst_weight $weight
-constraints:cst_fa_file input_files/$base.cst
-constraints:cst_fa_weight $weight

-out:file:fullatom
-out:pdb
-out:path:all output_files
-out:suffix _cst_relax
-overwrite
EOF
} > "$out/input_files/relax.options"

# Rosetta is needed to run, not to prepare and check. Everything above works
# without it, so a restraint file can be validated on a laptop.
require_rosetta relax
exe=$(rosetta_exe relax)
mpi_exe=$(rosetta_exe relax mpi 2>/dev/null || true)
if [[ -n "$mpi_exe" && "${CORES:-1}" -gt 1 ]]; then
    runline="mpirun -np ${CORES} $mpi_exe @input_files/relax.options"
else
    runline="$exe @input_files/relax.options"
    [[ "${CORES:-1}" -gt 1 ]] && warn "no MPI build configured; $nstruct structures will run on one core.
The original script requested 24 MPI ranks and then called the serial binary,
so 23 cores sat idle for the whole job. Set ROSETTA_MPI_SUFFIX to avoid that."
fi

cat > "$out/run.sh" <<EOF
#!/usr/bin/env bash
set -euo pipefail
${ROSETTA_SETUP:-}
cd "\$(dirname "\$0")"
$runline

shopt -s nullglob
sc=(output_files/*.sc output_files/*.fasc)
((\${#sc[@]})) || { echo "no score file written" >&2; exit 1; }

# Rank on total_score, then look at the restraint term separately: a structure
# can score well overall while leaving the restraints unsatisfied.
"\$QDM_ROOT/15-rosetta-modelling/extract/rosetta-scores.sh" \\
    -o output_files/scores -k total_score output_files
echo
echo "Restraint satisfaction, atom_pair_constraint by model:"
awk -F, 'NR==1{for(i=1;i<=NF;i++) if(\$i=="atom_pair_constraint") c=i; next}
         c{printf "  %-45s %10s\\n", \$1, \$c}' \\
    output_files/scores/scores-sorted.csv | head -10
EOF
chmod +x "$out/run.sh"

info "wrote $out/input_files/relax.options (weight $weight, nstruct $nstruct)"
if ((queue)); then
    printf '\nSubmit it with:\n  bin/qdm-submit %s\n' "$out"
else
    QDM_ROOT="$QDM_ROOT" "$out/run.sh"
fi
