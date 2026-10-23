#!/usr/bin/env bash
# =============================================================================
#  orca-input.sh - ORCA 5/6 jobs for a set of molecules, one folder each
#
#  Each step is a separate input (NAME-STEP.inp) so every output stays small and
#  easy to parse; run.sh runs the steps in order. Steps after "opt" start from
#  the optimised geometry (NAME-opt.xyz) and read its orbitals (MORead).
#  def2 basis sets carry their own ECPs (Rb onwards), so metal complexes need
#  no extra basis input.
#
#  usage: orca-input.sh [options] [-l list -d dir | files...]
#    -j STEPS  comma separated (default opt):
#       opt      geometry optimisation + frequencies (thermochemistry, IR)
#       sp       single point             freq   frequencies only
#       temps    thermochemistry at several temperatures            [-T]
#       td       TD-DFT singlets and triplets (UV-Vis)               [-t]
#       soc      TD-DFT with spin-orbit coupling (closed shell only) [-t]
#       nmr      NMR shieldings
#       excited  S0 and T1 optimisations + Hessians (NAME-GS, NAME-TG), the
#                starting point for fluorescence/phosphorescence rates with
#                orca-esd-input.sh
#       casscf   three steps: SCF orbitals -> CASSCF -> NEVPT2 (+ SOC with -S)
#                                                           [-a nel,norb -r roots -u mults]
#       neb      NEB-TS between NAME-r and NAME-p (+ Freq on the TS)  [-i images]
#    -m METHOD   functional / method  (default B3LYP)
#    -D DISP     D3BJ, D4 or none (default D3BJ)
#    -b BASIS    (default def2-TZVP); RIJCOSX with def2/J is added for hybrids
#    -s SOLV     smd:SOLVENT or cpcm:SOLVENT, e.g. smd:methanol (default gas)
#    -T LIST     temperatures for "temps" (default "268.15 273.15 283.15 293.15 298.15 303.15 308.15 313.15 318.15")
#    -t N        TD-DFT roots (default 25)
#    -a NEL,NORB active space for casscf (default 6,5)
#    -r ROOTS    CASSCF roots per multiplicity (default 2)
#    -u MULTS    multiplicities for the NEVPT2 step, e.g. "1,3" with -r "4,3"
#    -S          include spin-orbit coupling (QDPT) in the NEVPT2 step
#    -i N        NEB images (default 8)
#    -c FILE     charge/multiplicity table (default chargemult.txt)
#    -n NPROC -M MEM (total, e.g. 32G) -o OUTDIR (default orca-jobs)
#  example:
#    orca-input.sh -l list -d xyz -j opt,td,nmr -s smd:methanol
#    orca-input.sh -l reactions -d xyz -j neb          # needs NAME-r and NAME-p
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"
source "$QDM_ROOT/lib/jobscript.sh"
load_cluster_conf
steps=opt; method=B3LYP; disp=D3BJ; basis=def2-TZVP; solv=""; nroots=25
temps="268.15 273.15 283.15 293.15 298.15 303.15 308.15 313.15 318.15"
active=6,5; roots=2; mults=""; socnev=0; nimages=8; cmfile=chargemult.txt
nproc=${CORES:-8}; mem=${MEM:-16G}; outdir=orca-jobs; list=""; indir=.
while getopts "j:m:D:b:s:T:t:a:r:u:Si:c:n:M:o:l:d:h" o; do
    case $o in
        j) steps=$OPTARG ;; m) method=$OPTARG ;; D) disp=$OPTARG ;; b) basis=$OPTARG ;;
        s) solv=$OPTARG ;; T) temps=$OPTARG ;; t) nroots=$OPTARG ;; a) active=$OPTARG ;;
        r) roots=$OPTARG ;; u) mults=$OPTARG ;; S) socnev=1 ;; i) nimages=$OPTARG ;;
        c) cmfile=$OPTARG ;; n) nproc=$OPTARG ;; M) mem=$OPTARG ;; o) outdir=$OPTARG ;;
        l) list=$OPTARG ;; d) indir=$OPTARG ;; *) sed -n '2,45p' "$0"; exit 0 ;;
    esac
done
shift $((OPTIND-1))
CORES=$nproc; MEM=$mem; JOB_THREADED=0         # ORCA runs MPI processes
# %maxcore is memory per process in MB; keep 25 % headroom as ORCA recommends
memmb=$(echo "$mem" | awk '{v=$1+0; u=toupper($1); if (u ~ /T/) v*=1024*1024; else if (u ~ /G/) v*=1024; print v}')
maxcore=$(( memmb * 3 / 4 / nproc ))

# ---------- keyword pieces ----------
dk=""; [[ "${disp,,}" != none ]] && dk=" $disp"
ri=""
case "${method^^}" in
    B3LYP|PBE0|M06|M062X|TPSSH|CAM-B3LYP|WB97X*|B2PLYP*|PW6B95|BHANDHLYP) ri=" RIJCOSX def2/J" ;;
    BP86|PBE|BLYP|TPSS|R2SCAN|B97-D3|REVPBE) ri=" RI def2/J" ;;
esac
solkw=""; solblock=""
if [[ -n "$solv" ]]; then
    model=${solv%%:*}; sname=${solv#*:}
    solkw=" CPCM($sname)"
    case "${model,,}" in
        smd)  solblock=$'%cpcm\n  smd true\n  SMDsolvent "'"$sname"$'"\nend' ;;
        cpcm) ;;
        *) die "solvent model must be smd or cpcm" ;;
    esac
fi
header() {   # header "extra ! keywords" -> common input head
    echo "! $method$dk $basis$ri TightSCF$solkw $1"
    echo "%pal nprocs $nproc end"
    echo "%maxcore $maxcore"
    [[ -n "$solblock" ]] && echo "$solblock"
    echo "%scf MaxIter 500 end"
    return 0
}
geom() { echo "* xyzfile $1 $2 $3"; }     # geom CHARGE MULT FILE
# CASSCF steps start from Hartree-Fock orbitals and use RI-JK (no functional)
header_cas() {
    echo "! $basis def2/JK RIJK TightSCF$solkw $1"
    echo "%pal nprocs $nproc end"
    echo "%maxcore $maxcore"
    [[ -n "$solblock" ]] && echo "$solblock"
    echo "%scf MaxIter 500 end"
    return 0
}

tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
mkdir -p "$outdir"
mapfile -t names < <( if [[ -n "$list" ]]; then read_list "$list"; fi; for f in "$@"; do strip_ext "$f"; done )
((${#names[@]})) || die "no molecules (use -l list -d dir or file names)"
for name in "${names[@]}"; do
    d="$outdir/$name"; mkdir -p "$d"
    q=$(get_charge "$name" "$cmfile"); m=$(get_mult "$name" "$cmfile")
    cmds=(); src=""
    if [[ ",$steps," != *",neb,"* ]]; then
        src=$(find_structure "$name" "$indir") || { warn "no structure for $name"; continue; }
        to_format "$src" "$d/$name.xyz"
        check_multiplicity "$d/$name.xyz" "$q" "$m"
    fi
    g="$name.xyz"; mo=""                     # geometry and orbitals for the next step
    for s in ${steps//,/ }; do
        inp="$d/$name-$s.inp"
        case $s in
        opt)   { header "Opt Freq"; echo "%geom MaxIter 500 end"; geom "$q" "$m" "$g"; } > "$inp"
               g="$name-opt.xyz"; mo="$name-opt.gbw" ;;
        sp)    { header ""; geom "$q" "$m" "$g"; } > "$inp" ;;
        freq)  { header "Freq"; geom "$q" "$m" "$g"; } > "$inp" ;;
        temps) { header "Freq${mo:+ MORead}"; [[ -n "$mo" ]] && echo "%moinp \"$mo\""
                 echo "%freq Temp $(echo $temps | sed 's/ /, /g') end"; geom "$q" "$m" "$g"; } > "$inp" ;;
        td)    { header "${mo:+MORead}"; [[ -n "$mo" ]] && echo "%moinp \"$mo\""
                 printf '%%tddft\n  NRoots %s\n  Triplets true\nend\n' "$nroots"; geom "$q" "$m" "$g"; } > "$inp" ;;
        soc)   [[ $m == 1 ]] || { warn "$name: SOC-TDDFT needs a closed shell singlet, step skipped"; continue; }
               { header "RI-SOMF(1X)${mo:+ MORead}"; [[ -n "$mo" ]] && echo "%moinp \"$mo\""
                 printf '%%tddft\n  NRoots %s\n  DoSOC true\nend\n' "$nroots"; geom "$q" "$m" "$g"; } > "$inp" ;;
        nmr)   { header "NMR${mo:+ MORead}"; [[ -n "$mo" ]] && echo "%moinp \"$mo\""; geom "$q" "$m" "$g"; } > "$inp" ;;
        excited)
               [[ $m == 1 ]] || { warn "$name: 'excited' expects a singlet ground state"; continue; }
               tm=3
               { header "Opt Freq"; echo "%geom MaxIter 500 end"; geom "$q" "$m" "$g"; } > "$d/$name-GS.inp"
               { header "Opt Freq"; echo "%geom MaxIter 500 end"; geom "$q" "$tm" "$g"; } > "$d/$name-TG.inp"
               cmds+=("\$ORCA_EXE $name-GS.inp > $name-GS.out" "\$ORCA_EXE $name-TG.inp > $name-TG.out")
               continue ;;
        casscf)
               nel=${active%,*}; norb=${active#*,}; um=${mults:-$m}; ur=${roots}
               { header_cas "NormalPrint"; echo "%base \"$name-cas1\""; geom "$q" "$m" "$g"; } > "$d/$name-cas1.inp"
               { header_cas "MORead"; echo "%moinp \"$name-cas1.gbw\""; echo "%base \"$name-cas2\""
                 printf '%%casscf\n  nel %s\n  norb %s\n  mult %s\n  nroots %s\nend\n' "$nel" "$norb" "$m" "${roots%%,*}"
                 geom "$q" "$m" "$g"; } > "$d/$name-cas2.inp"
               { header_cas "MORead"; echo "%moinp \"$name-cas2.gbw\""; echo "%base \"$name-cas3\""
                 printf '%%casscf\n  nel %s\n  norb %s\n  mult %s\n  nroots %s\n  PTMethod SC_NEVPT2\n' "$nel" "$norb" "$um" "$ur"
                 ((socnev)) && printf '  rel\n    dosoc true\n  end\n'
                 echo "end"; geom "$q" "$m" "$g"; } > "$d/$name-cas3.inp"
               for k in 1 2 3; do cmds+=("\$ORCA_EXE $name-cas$k.inp > $name-cas$k.out"); done
               continue ;;
        neb)
               r=$(find_structure "$name-r" "$indir") && p=$(find_structure "$name-p" "$indir") ||
                   { warn "$name: NEB needs $name-r and $name-p"; continue; }
               to_format "$r" "$d/$name-r.xyz"; to_format "$p" "$d/$name-p.xyz"
               er=$(awk 'NR>2{printf $1}' "$d/$name-r.xyz"); ep=$(awk 'NR>2{printf $1}' "$d/$name-p.xyz")
               [[ "$er" == "$ep" ]] || { warn "$name: reactant and product atoms are not in the same order"; continue; }
               { header "NEB-TS Freq"
                 printf '%%neb\n  NEB_End_XYZFile "%s"\n  NImages %s\n  PreOpt_Ends false\nend\n' "$name-p.xyz" "$nimages"
                 geom "$q" "$m" "$name-r.xyz"; } > "$inp" ;;
        *) die "unknown step: $s" ;;
        esac
        cmds+=("\$ORCA_EXE $name-$s.inp > $name-$s.out")
    done
    # ORCA must be called with its full path for parallel runs (ORCA_EXE in cluster.conf)
    write_jobscript "o-$name" "$PWD/$d" ORCA_SETUP "ORCA_EXE=\"${ORCA_EXE:-orca}\"" "${cmds[@]}" > "$d/run.sh"
    chmod +x "$d/run.sh"
    info "$name: q=$q mult=$m steps=$steps -> $d/"
done
