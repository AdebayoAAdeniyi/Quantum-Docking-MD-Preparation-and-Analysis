#!/usr/bin/env bash
# =============================================================================
#  common.sh - shared helpers for every workflow in this repository
#
#  Source it at the top of a script:
#      source "$(dirname "$0")/../../lib/common.sh"
#
#  Provides: logging, dependency checks, molecule lists, structure conversion
#  through Open Babel, element handling, template rendering, unit constants
#  and a single awk matrix transpose (the old scripts each wrote their own).
# =============================================================================

# Repository root (works no matter where the calling script lives)
QDM_ROOT="${QDM_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
export QDM_ROOT

# ----------------------------------------------------------------------------
# Unit conversions (CODATA 2018)
# ----------------------------------------------------------------------------
HARTREE2EV=27.211386245988
HARTREE2KJ=2625.4996394799
HARTREE2KCAL=627.50947406
HARTREE2CM=219474.6313632
EV2KJ=96.48533212
KB_HARTREE=3.166811563e-6        # Boltzmann constant in Hartree/K
R_KJ=0.008314462618              # gas constant kJ/(mol K)

# ----------------------------------------------------------------------------
# Logging
# ----------------------------------------------------------------------------
if [[ -t 2 ]]; then _c_g=$'\e[32m'; _c_y=$'\e[33m'; _c_r=$'\e[31m'; _c_0=$'\e[0m'
else _c_g=""; _c_y=""; _c_r=""; _c_0=""; fi
info()  { printf '%s[info]%s %s\n'  "$_c_g" "$_c_0" "$*" >&2; }
warn()  { printf '%s[warn]%s %s\n'  "$_c_y" "$_c_0" "$*" >&2; }
die()   { printf '%s[error]%s %s\n' "$_c_r" "$_c_0" "$*" >&2; exit 1; }

# require cmd1 cmd2 ... : stop with a clear message when a program is missing
require() {
    local c missing=()
    for c in "$@"; do command -v "$c" >/dev/null 2>&1 || missing+=("$c"); done
    ((${#missing[@]} == 0)) || die "missing program(s): ${missing[*]} (see install/README.md)"
}

# ----------------------------------------------------------------------------
# Molecule lists
#   read_list FILE      -> one name per line, comments (#) and blanks removed,
#                          extensions stripped, so 'aspirin.xyz' -> 'aspirin'
#   gather_structures   -> structure files from "-l list -d dir" and/or file names
# ----------------------------------------------------------------------------
read_list() {
    [[ -f "$1" ]] || die "list file not found: $1"
    sed -e 's/#.*//' -e 's/[[:space:]]*$//' -e '/^[[:space:]]*$/d' "$1" |
        awk '{print $1}'
}

# gather_structures LIST DIR [files...]: print structure paths, one per line,
# from a list file (names looked up in DIR) and/or explicit file names
gather_structures() {
    local list="$1" dir="$2" n f; shift 2
    for f in "$@"; do [[ -f "$f" ]] && echo "$f" || warn "not found: $f"; done
    [[ -z "$list" ]] && return 0
    while read -r n; do
        f=$(find_structure "$n" "$dir") && echo "$f" || warn "no structure file for '$n' in $dir"
    done < <(read_list "$list")
}

strip_ext() { local b; b=$(basename "$1"); echo "${b%.*}"; }
ext_of()    { local b; b=$(basename "$1"); echo "${b##*.}" | tr 'A-Z' 'a-z'; }

# find_structure NAME [DIR]: locate NAME.{xyz,pdb,mol2,sdf,mol,cml,smi,gjf,com,log,out}
find_structure() {
    local name="$1" dir="${2:-.}" e
    [[ -f "$name" ]] && { echo "$name"; return 0; }
    for e in xyz pdb mol2 sdf mol cml smi gjf com log out; do
        [[ -f "$dir/$name.$e" ]] && { echo "$dir/$name.$e"; return 0; }
    done
    return 1
}

# ----------------------------------------------------------------------------
# Structure conversion with Open Babel
#   to_format IN OUT [extra obabel options]
#   SMILES or 2D input gets 3D coordinates and a short MMFF94 clean-up.
# ----------------------------------------------------------------------------
to_format() {
    local in="$1" out="$2"; shift 2
    require obabel
    local ie oe; ie=$(ext_of "$in"); oe=$(ext_of "$out")
    case "$ie" in log|out) ie=g09 ;; com|gjf) ie=gjf ;; esac
    local gen=()
    if [[ "$ie" == smi ]] || { [[ "$ie" =~ ^(sdf|mol)$ ]] && is_2d "$in"; }; then
        gen=(--gen3d --minimize --ff MMFF94 --steps 500)
    fi
    obabel -i"$ie" "$in" -o"$oe" -O "$out" "${gen[@]}" "$@" 2>/dev/null ||
        die "Open Babel could not convert $in -> $out"
    [[ -s "$out" ]] || die "conversion produced an empty file: $out"
}

# is_2d FILE.sdf: true when every z coordinate in the first molecule is zero
is_2d() {
    awk 'NR==4{n=$1} NR>4 && NR<=4+n {if ($3+0 != 0) nz=1} END{exit nz?1:0}' "$1"
}

# xyz_body FILE: "El x y z" lines for the first structure in any format
xyz_body() {
    local tmp; tmp=$(mktemp --suffix=.xyz)
    to_format "$1" "$tmp"
    awk 'NR==1{n=$1} NR>2 && NR<=n+2 {printf "%-3s %14.8f %14.8f %14.8f\n",$1,$2,$3,$4}' "$tmp"
    rm -f "$tmp"
}

# elements FILE: unique element symbols, in order of first appearance
elements() { xyz_body "$1" | awk '!seen[$1]++{print $1}'; }

# Transition metals, lanthanides and heavy main-group atoms that usually need
# an effective core potential (ECP) basis set such as LanL2DZ, SDD or def2-ECP.
QDM_HEAVY="Sc Ti V Cr Mn Fe Co Ni Cu Zn Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe
Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi"
is_heavy() { [[ " ${QDM_HEAVY//$'\n'/ } " == *" $1 "* ]]; }

# split_elements FILE -> sets LIGHT_ELEMENTS and HEAVY_ELEMENTS (space separated)
split_elements() {
    LIGHT_ELEMENTS=""; HEAVY_ELEMENTS=""
    local e
    for e in $(elements "$1"); do
        if is_heavy "$e"; then HEAVY_ELEMENTS+="$e "; else LIGHT_ELEMENTS+="$e "; fi
    done
    LIGHT_ELEMENTS="${LIGHT_ELEMENTS% }"; HEAVY_ELEMENTS="${HEAVY_ELEMENTS% }"
}

# metals_first FILE: coordinates with heavy atoms first and hydrogens last.
# (The old prepare-new scripts did this with one grep per element.)
metals_first() {
    xyz_body "$1" | awk -v heavy="$QDM_HEAVY" '
        BEGIN{n=split(heavy,h," "); for(i=1;i<=n;i++) H[h[i]]=1}
        { if ($1 in H) a[++na]=$0; else if ($1=="H") c[++nc]=$0; else b[++nb]=$0 }
        END{ for(i=1;i<=na;i++)print a[i]; for(i=1;i<=nb;i++)print b[i]; for(i=1;i<=nc;i++)print c[i] }'
}

# natoms FILE
natoms() { xyz_body "$1" | wc -l; }

# total_electrons FILE CHARGE: used to check a charge/multiplicity pair
total_electrons() {
    local z; z=$(xyz_body "$1" | awk '{print $1}' | awk -f "$QDM_ROOT/lib/awk/atomic_number.awk" | awk '{s+=$1}END{print s}')
    echo $(( z - ${2:-0} ))
}

# check_multiplicity FILE CHARGE MULT: warn when the pair is impossible
check_multiplicity() {
    local ne; ne=$(total_electrons "$1" "$2")
    if (( (ne + $3) % 2 == 0 )); then
        die "$(basename "$1"): charge $2 gives $ne electrons, multiplicity $3 is impossible"
    fi
}

# ----------------------------------------------------------------------------
# Charge / multiplicity table
#   A plain text file with "name charge multiplicity" per line, e.g.
#       aspirin   0 1
#       Fe-hfaa3  0 6
#   get_charge NAME [table]   get_mult NAME [table]   (defaults 0 and 1)
# ----------------------------------------------------------------------------
get_charge() { awk -v n="$1" '$1==n{print $2; f=1} END{if(!f) print 0}' "${2:-chargemult.txt}" 2>/dev/null || echo 0; }
get_mult()   { awk -v n="$1" '$1==n{print $3; f=1} END{if(!f) print 1}' "${2:-chargemult.txt}" 2>/dev/null || echo 1; }

# ----------------------------------------------------------------------------
# Templates: render FILE replacing {{KEY}} with the value of variable KEY
#   render_template file KEY1 KEY2 ...   (values taken from the environment)
# ----------------------------------------------------------------------------
render_template() {
    local tpl="$1"; shift
    [[ -f "$tpl" ]] || die "template not found: $tpl"
    (   # export the keys only inside this subshell so awk can read them
        for k in "$@"; do export "$k"; done
        awk -v keys="$*" '
            # literal (not regex) replacement, safe for "&", "\" and "/"
            function repl(s, pat, val,   out, i) {
                out = ""
                while ((i = index(s, pat)) > 0) { out = out substr(s, 1, i-1) val; s = substr(s, i+length(pat)) }
                return out s
            }
            BEGIN{n=split(keys,K," ")}
            { line=$0; for(i=1;i<=n;i++) line=repl(line, "{{" K[i] "}}", ENVIRON[K[i]]); print line }' "$tpl"
    )
}

# ----------------------------------------------------------------------------
# Small numeric helpers (awk based, no bc needed)
# ----------------------------------------------------------------------------
calc() { awk "BEGIN{printf \"%.10g\n\", $*}"; }

# transpose: rows <-> columns of a whitespace table (one copy for the repo)
transpose() { awk -f "$QDM_ROOT/lib/awk/transpose.awk" "$@"; }

# csv_join a b c -> "a,b,c"
csv_join() { local IFS=,; echo "$*"; }

# last_match PATTERN FILE: print the last line that matches
last_match() { grep -E "$1" "$2" | tail -n 1; }

# ----------------------------------------------------------------------------
# Rosetta
#   rosetta_exe APP [mpi]  -> full command for a Rosetta application
#
# Rosetta binaries are named APP.BUILD.OSCOMPILERMODE, for example
# relax.default.linuxgccrelease. The build part varies between installations,
# so it comes from config/cluster.conf rather than being written into scripts.
# ----------------------------------------------------------------------------
rosetta_exe() {
    local app="$1" kind="${2:-serial}" suffix exe
    if [[ "$kind" == mpi ]]; then suffix="${ROSETTA_MPI_SUFFIX:-mpi.linuxgccrelease}"
    else suffix="${ROSETTA_SUFFIX:-default.linuxgccrelease}"; fi
    for exe in "${ROSETTA_BIN:+$ROSETTA_BIN/}$app.$suffix" "$app.$suffix" "$app"; do
        command -v "$exe" >/dev/null 2>&1 && { echo "$exe"; return 0; }
        [[ -x "$exe" ]] && { echo "$exe"; return 0; }
    done
    return 1
}

# require_rosetta APP ... : stop with a message naming the missing application
require_rosetta() {
    local app missing=()
    for app in "$@"; do rosetta_exe "$app" >/dev/null || missing+=("$app"); done
    ((${#missing[@]} == 0)) || die "Rosetta application(s) not found: ${missing[*]}
Set ROSETTA_BIN and ROSETTA_SUFFIX in config/cluster.local.conf.
Rosetta is free for non-commercial users; see docs/third-party-material.md."
}
