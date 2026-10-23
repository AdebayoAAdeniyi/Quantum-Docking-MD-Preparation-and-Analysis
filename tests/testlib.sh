#!/usr/bin/env bash
# Small helpers shared by the test scripts.
#   fetch URL DEST           download once into tests/data (cached)
#   check "label" ACTUAL EXPECTED [tolerance]   numeric or string comparison
#   need cmd                  skip the test file when a program is missing
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
finish() { echo "  $PASS passed, $FAIL failed"; ((FAIL == 0)); }
