#!/usr/bin/env bash
# Small helpers shared by the test scripts.
#   fetch URL DEST           download once into tests/data (cached)
#   check "label" ACTUAL EXPECTED [tolerance]   numeric or string comparison
#   need cmd                  skip the test file when a program is missing
#   fixture NAME              path to tests/fixtures/NAME, or "" with a warning
TESTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
QDM_ROOT="$(cd "$TESTS_DIR/.." && pwd)"; export QDM_ROOT
DATA="$TESTS_DIR/data"; mkdir -p "$DATA"
PASS=0; FAIL=0
CCLIB=https://raw.githubusercontent.com/cclib/cclib/master/data
CCLIBDATA=https://raw.githubusercontent.com/cclib/cclib-data/master
fetch() {
    local url=$1 dest="$DATA/$2"
    [[ -s "$dest" ]] && return 0
    mkdir -p "$(dirname "$dest")"
    curl -fsSL "$url" -o "$dest" || { echo "  cannot download $url"; rm -f "$dest"; return 1; }
}
check() {
    local label=$1 got=$2 want=$3 tol=${4:-}
    if [[ -n "$tol" ]] && awk -v a="$got" -v b="$want" -v t="$tol" 'BEGIN{d=a-b; if(d<0)d=-d; exit !(a!="" && d<=t)}'; then ok=1
    elif [[ -z "$tol" && "$got" == "$want" ]]; then ok=1; else ok=0; fi
    if ((ok)); then PASS=$((PASS+1)); printf '  ok    %s\n' "$label"
    else FAIL=$((FAIL+1)); printf '  FAIL  %s: got "%s", expected "%s"\n' "$label" "$got" "$want"; fi
}
need() { command -v "$1" >/dev/null || { echo "  skipped: $1 not installed"; exit 0; }; }
# A fixture is a fragment of real program output committed for the parser tests.
# If one is absent, say so plainly and let the rest of the file run: a missing
# fixture is a packaging problem, not a parser failure, and an awk crash hides
# the difference.
fixture() {
    local f="$TESTS_DIR/fixtures/$1"
    if [[ -s "$f" ]]; then printf '%s' "$f"; return 0; fi
    printf '  missing fixture tests/fixtures/%s, checks that need it skipped\n' "$1" >&2
    return 1
}
finish() { echo "  $PASS passed, $FAIL failed"; ((FAIL == 0)); }
