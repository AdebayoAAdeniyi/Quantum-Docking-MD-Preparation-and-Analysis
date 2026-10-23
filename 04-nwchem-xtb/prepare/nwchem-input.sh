#!/usr/bin/env bash
# =============================================================================
#  nwchem-input.sh - NWChem DFT jobs for a set of molecules or reactions
#
#  usage: nwchem-input.sh -j JOB [options] [-l list -d dir | files...]
#    -j opt    optimisation + numerical frequencies, thermochemistry at -T temperatures
#       sp     single point energy
#       ts     transition state search (task dft saddle) from NAME-ts, then frequencies
#       neb    nudged elastic band between NAME-r and NAME-p (same atom order);
#              NAME-ts, if present, is used as the middle of the starting path (xyz_path)
#       bsse   counterpoise corrected interaction energy; fragments with -f "1-10/11-14"
#              (or per molecule in the list, second column)
#    -x XC     functional (default b3lyp)      -b BASIS   (default 6-311+g(2d,p))
#    -D        add D3(BJ) dispersion           -s SOLV    COSMO solvent name or dielectric
#    -T LIST   temperatures in K (default "200 230 250 270 298.15 300 320 350 380 400")
#    -B N      NEB beads (default 20)          -c FILE    charge/multiplicity table
#    -n NPROC  -M MEM  -o OUTDIR (default nwchem-jobs)
#  The original scripts kept the molecule names and charges inside each script;
#  here they come from a list and chargemult.txt.
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"
source "$QDM_ROOT/lib/jobscript.sh"
load_cluster_conf
job=""; xc=b3lyp; basis="6-311+g(2d,p)"; disp=0; solv=""; temps="200 230 250 270 298.15 300 320 350 380 400"
nbeads=20; frags=""; cmfile=chargemult.txt; nproc=${CORES:-8}; mem=${MEM:-16G}; outdir=nwchem-jobs; list=""; indir=.
while getopts "j:x:b:Ds:T:B:f:c:n:M:o:l:d:h" o; do
    case $o in j) job=$OPTARG ;; x) xc=$OPTARG ;; b) basis=$OPTARG ;; D) disp=1 ;; s) solv=$OPTARG ;; T) temps=$OPTARG ;;
               B) nbeads=$OPTARG ;; f) frags=$OPTARG ;; c) cmfile=$OPTARG ;; n) nproc=$OPTARG ;; M) mem=$OPTARG ;;
               o) outdir=$OPTARG ;; l) list=$OPTARG ;; d) indir=$OPTARG ;; *) sed -n '2,21p' "$0"; exit 0 ;; esac
done
shift $((OPTIND-1))
[[ "$job" =~ ^(opt|sp|ts|neb|bsse)$ ]] || die "choose -j opt, sp, ts, neb or bsse"
CORES=$nproc; MEM=$mem
memmb=$(echo "$mem" | awk '{v=$1+0; u=toupper($1); if (u ~ /T/) v*=1048576; else if (u ~ /G/) v*=1024; print v}')
permb=$(( memmb * 3 / 4 / nproc ))                     # NWChem "memory" is per process
expand() { echo "$1" | tr ',' '\n' | awk -F- 'NF==2{for(i=$1;i<=$2;i++) printf "%d ", i; next} {printf "%d ", $1}'; }

head_nw() {   # head_nw NAME CHARGE
    cat <<EOF
title "$1 $job"
echo
start $1
memory total $permb mb
charge $2
permanent_dir ./perm
scratch_dir ./perm
EOF
}
dft_block() {   # dft_block MULT
    echo "basis spherical"; echo "  * library $basis"
    # counterpoise: ghost atoms (bqX) need the same functions as the real atoms
    if [[ $job == bsse ]]; then for e in $GHOSTS; do echo "  bq$e library $e $basis"; done; fi
    echo "end"
    echo "dft"; echo "  xc $xc"; echo "  mult $1"; echo "  grid fine"; echo "  maxiter 500"
    echo "  convergence energy 1d-8"; ((disp)) && echo "  disp vdw 4"; echo "end"
    [[ -n "$solv" ]] && { if [[ "$solv" =~ ^[0-9.]+$ ]]; then echo "cosmo; dielec $solv; end"; else echo "cosmo; solvent $solv; end"; fi; }
    return 0
}
geom_block() {  # geom_block FILE [keyword]
    echo "geometry ${2:-} units angstrom nocenter noautosym noautoz"
    xyz_body "$1" | awk -v b="${3:-}" '{ printf "  %-14s %14.8f %14.8f %14.8f\n", $1 b, $2, $3, $4 }'
    echo "end"
}
freq_block() {
    local t=($temps)
    echo "freq"; echo "  temp ${#t[@]} ${temps}"; echo "end"
}

mkdir -p "$outdir"
while read -r name rest; do
    q=$(get_charge "$name" "$cmfile"); m=$(get_mult "$name" "$cmfile")
    d="$outdir/$name"; mkdir -p "$d/perm"; inp="$d/$name-$job.nw"
    case $job in
    opt|sp)
        f=$(find_structure "$name" "$indir") || { warn "no structure for $name"; continue; }
        check_multiplicity "$f" "$q" "$m"
        { head_nw "$name-$job" "$q"; geom_block "$f"; dft_block "$m"
          if [[ $job == opt ]]; then
              printf 'driver\n  maxiter 300\n  tight\n  xyz %s-opt\nend\n' "$name"
              echo "task dft optimize"; freq_block; echo "task dft freq"
          else echo "task dft energy"; fi; } > "$inp" ;;
    ts)
        f=$(find_structure "$name-ts" "$indir") || { warn "no $name-ts structure"; continue; }
        { head_nw "$name-ts" "$q"; geom_block "$f"; dft_block "$m"
          printf 'driver\n  maxiter 300\n  inhess 2\n  xyz %s-ts\nend\n' "$name"
          echo "task dft freq"                 # Hessian first, so the saddle search follows the right mode
          echo "task dft saddle"; freq_block; echo "task dft freq"; } > "$inp" ;;
    neb)
        r=$(find_structure "$name-r" "$indir") && p=$(find_structure "$name-p" "$indir") ||
            { warn "$name: needs $name-r and $name-p"; continue; }
        [[ "$(xyz_body "$r" | awk '{printf $1}')" == "$(xyz_body "$p" | awk '{printf $1}')" ]] ||
            { warn "$name: reactant and product atom order differ"; continue; }
        mid=$(find_structure "$name-ts" "$indir") || mid=""
        pathkw=""
        if [[ -n "$mid" ]]; then
            # a straight line between the end points can push atoms through each other; with a
            # TS guess the starting path is interpolated in two halves, R -> TS -> P
            { xyz_body "$r"; echo; xyz_body "$mid"; echo; xyz_body "$p"; } | awk -v nb="$nbeads" -v nm="$name" 'BEGIN { s = 0 }
                /^$/ { s++; next } { k = s; n[k]++; i = n[k]; e[i] = $1; X[k,i] = $2; Y[k,i] = $3; Z[k,i] = $4 }
                END { na = n[0]; h = int((nb - 1) / 2)
                      for (b = 0; b < nb; b++) {
                          if (b <= h) { a = 0; c = 1; t = b / h } else { a = 1; c = 2; t = (b - h) / (nb - 1 - h) }
                          print na; print nm " bead " b + 1
                          for (i = 1; i <= na; i++) printf "%-3s %14.8f %14.8f %14.8f\n", e[i], X[a,i] + t * (X[c,i] - X[a,i]), Y[a,i] + t * (Y[c,i] - Y[a,i]), Z[a,i] + t * (Z[c,i] - Z[a,i]) } }' > "$d/perm/$name-path.xyz"   # NWChem reads it from permanent_dir
            pathkw="  xyz_path $name-path.xyz"
        fi
        { head_nw "$name-neb" "$q"; geom_block "$r"; geom_block "$p" endgeom
          dft_block "$m"
          printf 'neb\n  nbeads %s\n  kbeads 1.0\n  maxiter 50\n  stepsize 1.0\n  print_shift 1\n%s\nend\n' "$nbeads" "$pathkw"
          echo "task dft neb ignore"; } > "$inp" ;;
    bsse)
        f=$(find_structure "$name" "$indir") || { warn "no structure for $name"; continue; }
        fr=${frags:-$(echo "$rest" | awk '{print $1}')}
        [[ "$fr" == */* ]] || { warn "$name: give two fragments, -f \"1-10/11-14\""; continue; }
        a=$(expand "${fr%%/*}"); b=$(expand "${fr#*/}"); GHOSTS=$(elements "$f" | xargs)
        { head_nw "$name-bsse" "$q"; geom_block "$f"; dft_block "$m"
          printf 'bsse\n  mon frag1 %s\n  mon frag2 %s\nend\n' "$a" "$b"
          echo "task dft energy"; } > "$inp" ;;
    esac
    write_jobscript "nw-$name" "$PWD/$d" NWCHEM_SETUP "mkdir -p perm" \
        "\${MPIRUN:-mpirun} -np \$NCORES \${NWCHEM_EXE:-nwchem} $name-$job.nw > $name-$job.log" > "$d/run-$job.sh"
    chmod +x "$d/run-$job.sh"
    info "$name ($job): $inp"
done < <( { [[ -n "$list" ]] && sed -e 's/#.*//' -e '/^[[:space:]]*$/d' "$list"; for a in "$@"; do strip_ext "$a"; done; } )
