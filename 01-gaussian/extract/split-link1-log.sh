#!/usr/bin/env bash
# =============================================================================
#  split-link1-log.sh - cut a multi-step (--Link1--) Gaussian log into one file
#                       per job step, named after what the step did
#
#  usage: split-link1-log.sh FILE.log [OUTDIR]
#     writes OUTDIR/NAME.01-opt-freq.log, OUTDIR/NAME.02-pop.log, ...
#     (OUTDIR defaults to NAME-steps next to the log) and prints the names.
#
#  A step ends at its "Normal termination" or "Error termination" line, so the
#  split does not depend on the version specific Link1 banner.
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"
log=${1:?usage: split-link1-log.sh FILE.log [OUTDIR]}
[[ -f "$log" ]] || die "no such file: $log"
name=$(strip_ext "$log")
out=${2:-$(dirname "$log")/$name-steps}
mkdir -p "$out"
rm -f "$out/$name".[0-9][0-9]-*.log
awk -v out="$out" -v name="$name" '
    function label(r,   l) {             # short label from the route section
        r = tolower(r); l = ""
        if (r ~ /qst[23]/)                 l = l "-qst"
        else if (r ~ /opt[^ ]*[(,=]ts/)    l = l "-ts"
        else if (r ~ /opt/)                l = l "-opt"
        if (r ~ /freq/)    l = l "-freq"
        if (r ~ /irc/)     l = l "-irc"
        if (r ~ /td *[=(]|tda|cis/) l = l "-td"
        if (r ~ /nmr/)     l = l "-nmr"
        if (r ~ /polar/)   l = l "-polar"
        if (r ~ /epr/)     l = l "-epr"
        if (r ~ /nbo/)     l = l "-nbo"
        if (r ~ /chelpg|mk|hly/) l = l "-esp"
        else if (r ~ /pop/) l = l "-pop"
        if (r ~ /wfx|wfn/) l = l "-wfx"
        if (r ~ /counterpoise/) l = l "-cp"
        if (l == "") l = "-sp"
        return substr(l, 2)
    }
    function flush() {
        if (n == 0) return
        step++; f = sprintf("%s/%s.%02d-%s.log", out, name, step, label(route))
        for (i = 1; i <= n; i++) print buf[i] > f
        close(f); print f; n = 0; route = ""; inroute = 0
    }
    { buf[++n] = $0 }
    /^ #/ && route == "" { inroute = 1 }
    inroute { if ($0 ~ /^ ---/) inroute = 0; else route = route $0 }
    /(Normal|Error) termination of Gaussian/ { flush() }
    END { if (n > 5) flush() }        # an unfinished last step
' "$log"
