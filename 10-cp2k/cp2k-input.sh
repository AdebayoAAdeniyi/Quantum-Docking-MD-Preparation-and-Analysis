#!/usr/bin/env bash
# =============================================================================
#  cp2k-input.sh - CP2K inputs that start from an Amber system (module 08)
#
#  usage: cp2k-input.sh -d AMBER_SYSTEM_DIR -j mm|qmmm [-c COORDS.rst7] [options]
#    mm     classical MD with the Amber force field (FIST, SPME), NPT_I by default
#    qmmm   QM/MM MD: the -q selection is treated with DFT (PBE-D3/DZVP-MOLOPT-SR,
#           GPW) and the rest with the Amber force field, electrostatic embedding
#           (Gaussian expansion of the electrostatic potential, GEEP)
#    -c     coordinates (default: the last Amber restart found: prod, equil, density, solv.inpcrd)
#    -q     QM atoms as an Amber residue name (e.g. LIG) or a list of atom indices "1-25,40"
#    -Q     charge of the QM region (default 0)
#    -n     MD steps (default 10000)   -t  time step fs (default 1.0 for mm, 0.5 for qmmm)
#    -T     temperature (default 300)  -E  ensemble NPT_I, NVT or NVE (default NPT_I for mm, NVT for qmmm)
#    -P     add a PLUMED file (plumed.dat) to the run (see module 12)
#  Writes NAME-mm.inp or NAME-qmmm.inp and run-cp2k.sh in the system folder.
# =============================================================================
source "$(dirname "$0")/../lib/common.sh"
source "$QDM_ROOT/lib/jobscript.sh"; load_cluster_conf
d=""; job=mm; crd=""; qsel=""; qq=0; steps=10000; dt=""; T=300; ens=""; plumed=0
while getopts "d:j:c:q:Q:n:t:T:E:Ph" o; do
    case $o in d) d=$OPTARG ;; j) job=$OPTARG ;; c) crd=$OPTARG ;; q) qsel=$OPTARG ;; Q) qq=$OPTARG ;; n) steps=$OPTARG ;;
               t) dt=$OPTARG ;; T) T=$OPTARG ;; E) ens=$OPTARG ;; P) plumed=1 ;; *) sed -n '2,19p' "$0"; exit 0 ;; esac
done
[[ -f "$d/solv.prmtop" ]] || die "no solv.prmtop in $d"
[[ "$job" =~ ^(mm|qmmm)$ ]] || die "-j mm or qmmm"
if [[ -z "$crd" ]]; then for c in prod.rst7 equil.rst7 density.rst7 solv.inpcrd; do [[ -f "$d/$c" ]] && { crd=$c; break; }; done; fi
[[ -f "$d/$crd" ]] || die "no coordinates"
name=$(basename "$d")
# CP2K guesses each atom's element from the coordinate file, and Amber force field
# atom types are not element symbols, so the coordinates are converted to a PDB with
# cpptraj, which writes the element column (and the box as CRYST1).
require cpptraj
pdb="$d/${crd%.*}-cp2k.pdb"
printf 'parm %s\ntrajin %s\ntrajout %s pdb\nrun\n' "$d/solv.prmtop" "$d/$crd" "$pdb" |
    cpptraj > /dev/null 2>&1 || die "cpptraj could not convert $crd to PDB"
[[ -s "$pdb" ]] || die "no coordinates written"
read -r a b c al be ga <<< "$(awk '/^CRYST1/{print $2, $3, $4, $5, $6, $7; exit}' "$pdb")"
[[ -n "$ga" ]] || die "$crd has no box information"
crd=$(basename "$pdb")
oct=$(awk -v x="$al" 'BEGIN{print (x > 100 && x < 120) ? 1 : 0}')
: "${dt:=$( [[ $job == mm ]] && echo 1.0 || echo 0.5 )}"; : "${ens:=$( [[ $job == mm ]] && echo NPT_I || echo NVT )}"
inp="$d/$name-$job.inp"
prmtop_field() {   # four-character fields of a prmtop flag, one per line
    awk -v f="%FLAG $1" '$0 ~ "^" f {g = 1; getline; next} /^%FLAG/ {g = 0}
        g { for (i = 1; i <= length($0); i += 4) { s = substr($0, i, 4); gsub(/ /, "", s); if (s != "") print s } }' "$d/solv.prmtop"
}
prmtop_numbers() { awk -v f="%FLAG $1" '$0 ~ "^" f {g = 1; getline; next} /^%FLAG/ {g = 0} g {for (i = 1; i <= NF; i++) print $i}' "$d/solv.prmtop"; }

qm_atoms() {   # 1-based atom indices of the QM region
    if [[ "$qsel" =~ ^[0-9,-]+$ ]]; then
        echo "$qsel" | tr ',' '\n' | awk -F- 'NF==2{for(i=$1;i<=$2;i++) print i; next}{print $1}'
    else
        paste -d' ' <(prmtop_field RESIDUE_LABEL) <(prmtop_numbers RESIDUE_POINTER) |
            awk -v r="$qsel" -v nat="$(prmtop_numbers POINTERS | head -1)" '
                { L[NR] = $1; P[NR] = $2; n = NR }
                END { for (k = 1; k <= n; k++) if (L[k] == r) { e = (k < n) ? P[k+1] - 1 : nat; for (j = P[k]; j <= e; j++) print j } }'
    fi
}

{
cat <<EOI
&GLOBAL
  PROJECT $name-$job
  RUN_TYPE MD
  PRINT_LEVEL LOW
&END GLOBAL
&FORCE_EVAL
  METHOD $( [[ $job == mm ]] && echo FIST || echo QMMM )
  STRESS_TENSOR $( [[ $ens == NPT_I && $job == mm ]] && echo ANALYTICAL || echo NONE )
  &MM
    &FORCEFIELD
      PARMTYPE AMBER
      PARM_FILE_NAME solv.prmtop
      EI_SCALE14 0.8333333
      VDW_SCALE14 0.5
      &SPLINE
        EMAX_SPLINE 1.0E14
        RCUT_NB [angstrom] 10.0
      &END SPLINE
    &END FORCEFIELD
    &POISSON
      &EWALD
        EWALD_TYPE SPME
        ALPHA [angstrom^-1] 0.35
        GMAX $(awk -v a="$a" 'BEGIN{g=int(a); if (g%2==0) g++; print g}')
      &END EWALD
    &END POISSON
  &END MM
EOI
if [[ $job == qmmm ]]; then
    [[ -n "$qsel" ]] || die "give the QM region with -q"
    mapfile -t qa < <(qm_atoms)
    # CP2K needs an explicit &LINK for every covalent bond cut by the QM/MM boundary
    whole=$(paste -d' ' <(prmtop_numbers RESIDUE_POINTER) <(seq 1 "$(prmtop_numbers RESIDUE_POINTER | wc -l)") |
        awk -v qm="${qa[*]:-}" -v nat="$(prmtop_numbers POINTERS | head -1)" '
            { P[NR] = $1; n = NR } END { m = split(qm, Q, " "); for (i = 1; i <= m; i++) S[Q[i] + 0] = 1
                for (k = 1; k <= n; k++) { e = (k < n) ? P[k+1] - 1 : nat; inq = out = 0
                    for (j = P[k]; j <= e; j++) if (j in S) inq++; else out++
                    if (inq && out) partial++ }
                print (partial ? "no" : "yes") }')
    [[ $whole == yes ]] || warn "the QM selection cuts through at least one residue; add &LINK sections to the &QMMM block (see docs/troubleshooting.md)"
((${#qa[@]})) || die "QM selection '$qsel' matched no atoms; residues in this system: $(prmtop_field RESIDUE_LABEL | sort -u | grep -vxE 'WAT|HOH|Na\+|Cl-|K\+' | paste -sd, -)"
    # element of each QM atom from the Amber atomic numbers
    mapfile -t z < <(awk '/%FLAG ATOMIC_NUMBER/{f=1; getline; next} /%FLAG/{f=0} f{for(i=1;i<=NF;i++) print $i}' "$d/solv.prmtop")
    declare -A byel=()
    sym=(X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca)
    for i in "${qa[@]}"; do e=${sym[${z[$((i-1))]}]}; byel[$e]+="$i "; done
    # the QM cell must enclose the QM region with room to spare: extent + 2 x 6 A
    qmcell=$(awk -v idx="${qa[*]}" 'BEGIN { n = split(idx, I, " "); for (i = 1; i <= n; i++) W[I[i]] = 1 }
        /^(ATOM|HETATM)/ { a = substr($0, 7, 5) + 0; if (!(a in W)) next
            x = substr($0, 31, 8) + 0; y = substr($0, 39, 8) + 0; z = substr($0, 47, 8) + 0
            if (!m++) { x0 = x1 = x; y0 = y1 = y; z0 = z1 = z }
            if (x < x0) x0 = x; if (x > x1) x1 = x
            if (y < y0) y0 = y; if (y > y1) y1 = y
            if (z < z0) z0 = z; if (z > z1) z1 = z }
        END { d = x1 - x0; if (y1 - y0 > d) d = y1 - y0; if (z1 - z0 > d) d = z1 - z0
              printf "%.1f %.1f %.1f", d + 12, d + 12, d + 12 }' "$pdb")      # cubic: the wavelet Poisson solver needs it
cat <<EOI
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT
    POTENTIAL_FILE_NAME GTH_POTENTIALS
    CHARGE $qq
    &MGRID
      CUTOFF 400
      COMMENSURATE
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-10
    &END QS
    &POISSON
      PERIODIC NONE
      POISSON_SOLVER WAVELET
    &END POISSON
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF 1.0E-6
      MAX_SCF 50
      &OT
        MINIMIZER DIIS
        PRECONDITIONER FULL_SINGLE_INVERSE
      &END OT
      &OUTER_SCF
        EPS_SCF 1.0E-6
        MAX_SCF 10
      &END OUTER_SCF
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
      &VDW_POTENTIAL
        POTENTIAL_TYPE PAIR_POTENTIAL
        &PAIR_POTENTIAL
          TYPE DFTD3
          PARAMETER_FILE_NAME dftd3.dat
          REFERENCE_FUNCTIONAL PBE
        &END PAIR_POTENTIAL
      &END VDW_POTENTIAL
    &END XC
  &END DFT
  &QMMM
    &CELL
      ABC $qmcell
      PERIODIC NONE
    &END CELL
    ECOUPL GAUSS
    USE_GEEP_LIB 12
    CENTER EVERY_STEP
EOI
    for e in "${!byel[@]}"; do echo "    &QM_KIND $e"; echo "      MM_INDEX ${byel[$e]}"; echo "    &END QM_KIND"; done
    echo "  &END QMMM"
fi
cat <<EOI
  &SUBSYS
    &CELL
      ABC [angstrom] $a $b $c
      ALPHA_BETA_GAMMA $al $be $ga
    &END CELL
    &TOPOLOGY
      CONN_FILE_FORMAT AMBER
      CONN_FILE_NAME solv.prmtop
      COORD_FILE_FORMAT PDB
      COORD_FILE_NAME $crd
    &END TOPOLOGY
EOI
# CP2K builds its atomic kinds from the Amber atom names and force field types, and
# neither is an element symbol ("nu", "hx", "C12"), so the element of every distinct
# name and type is stated explicitly. Without this CP2K stops with
# "Unknown element for KIND".
paste -d' ' <(prmtop_field ATOM_NAME) <(prmtop_field AMBER_ATOM_TYPE) <(prmtop_numbers ATOMIC_NUMBER) |
    awk -v qm="${qa[*]:-}" 'BEGIN { split("X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe", S, " ")
                                   n = split(qm, Q, " "); for (i = 1; i <= n; i++) isqm[Q[i] + 0] = 1 }
         { e = S[$3 + 1]; if (e == "") next
           # a kind that contains QM atoms also needs a basis set and a pseudopotential
           if (isqm[FNR]) qmel[e] = 1
           for (k = 1; k <= 2; k++) { if (seen[$k]++) { if (isqm[FNR]) needs[$k] = 1; continue }
               kinds[++nk] = $k; elem[$k] = e; if (isqm[FNR]) needs[$k] = 1 } }
         END { # CP2K may also create a kind named after the element itself, which then
               # needs the basis set as well
               for (e in qmel) if (!seen[e]++) { kinds[++nk] = e; elem[e] = e; needs[e] = 1 }
               for (i = 1; i <= nk; i++) { k = kinds[i]
                   print "    &KIND " k; print "      ELEMENT " elem[k]
                   if (needs[k] || (k in qmel)) { print "      BASIS_SET DZVP-MOLOPT-SR-GTH"; print "      POTENTIAL GTH-PBE" }
                   print "    &END KIND" } }'

cat <<EOI
  &END SUBSYS
&END FORCE_EVAL
&MOTION
  &MD
    ENSEMBLE $ens
    STEPS $steps
    TIMESTEP $dt
    TEMPERATURE $T
    &THERMOSTAT
      TYPE CSVR
      &CSVR
        TIMECON [fs] 100.0
      &END CSVR
    &END THERMOSTAT
EOI
[[ $ens == NPT_I ]] && printf '    &BAROSTAT\n      PRESSURE [bar] 1.0\n      TIMECON [fs] 1000.0\n    &END BAROSTAT\n'
echo "  &END MD"
((plumed)) && printf '  &FREE_ENERGY\n    METHOD METADYN\n    &METADYN\n      USE_PLUMED .TRUE.\n      PLUMED_INPUT_FILE plumed.dat\n    &END METADYN\n  &END FREE_ENERGY\n'
cat <<EOI
  &PRINT
    &TRAJECTORY
      FORMAT DCD
      &EACH
        MD 100
      &END EACH
    &END TRAJECTORY
    &RESTART
      &EACH
        MD 1000
      &END EACH
    &END RESTART
  &END PRINT
&END MOTION
EOI
} > "$inp"
((oct)) && [[ $job == qmmm ]] && warn "truncated octahedral box: check that the QM cell fits"
write_jobscript "cp2k-$name" "$PWD/$d" CP2K_SETUP "\${MPIRUN:-mpirun} -np \$NCORES \${CP2K_EXE:-cp2k.psmp} -i $(basename "$inp") -o ${name}-${job}.out" > "$d/run-cp2k-$job.sh"
chmod +x "$d/run-cp2k-$job.sh"
info "$inp (${ens}, $steps x $dt fs$( [[ $job == qmmm ]] && echo ", ${#qa[@]} QM atoms" ))"
