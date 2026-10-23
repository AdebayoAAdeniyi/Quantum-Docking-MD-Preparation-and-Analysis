#!/usr/bin/env bash
# =============================================================================
#  cv-process.sh - cyclic voltammetry data: ferrocene referencing, onset
#                  potentials and frontier orbital energies
#
#  Reads the text export of the potentiostat (two numeric columns, potential and
#  current; header lines and "Segment" lines are skipped wherever they are).
#
#  usage: cv-process.sh -F ferrocene.TXT [options] analyte1.TXT analyte2.TXT ...
#    -F FILE    ferrocene (internal or external standard) voltammogram
#    -w LO:HI   potential window of the Fc/Fc+ couple (same units as the file);
#               default: the whole scan
#    -u mV|V    potential units in the files (default mV, as exported by the
#               instrument used in the original work)
#    -c N       cycle to use: a number or "last" (default last); cycles are cut
#               where the scan returns to its starting direction
#    -d red|ox  process the reduction (cathodic) or the oxidation (anodic) wave
#    -o LO:HI   window (vs Fc+/Fc, V) that contains the rising part of that wave;
#               default: the whole sweep
#    -b LO:HI   baseline window (vs Fc+/Fc, V) before the wave (default: the first
#               15 % of the sweep in the direction of the wave)
#    -E EV      energy of Fc+/Fc below vacuum (default 4.8 eV)
#    -g EG      optical gap in eV (from uv-onset.sh) to get the other orbital
#    -O DIR     output folder (default cv-results)
#  Method: E1/2(Fc) = (Epa + Epc)/2; analyte potentials are shifted to Fc+/Fc;
#  the onset is the crossing of the baseline tangent with the tangent at the
#  steepest point of the rising current; E_LUMO = -(E_red,onset + 4.8) eV,
#  E_HOMO = -(E_ox,onset + 4.8) eV (and E_HOMO = E_LUMO - Eg_opt with -g).
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"
fc=""; fw=""; unit=mV; cyc=last; dir=red; ow=""; bw=""; evac=4.8; eg=""; out=cv-results
while getopts "F:w:u:c:d:o:b:E:g:O:h" o; do
    case $o in F) fc=$OPTARG ;; w) fw=$OPTARG ;; u) unit=$OPTARG ;; c) cyc=$OPTARG ;; d) dir=$OPTARG ;;
               o) ow=$OPTARG ;; b) bw=$OPTARG ;; E) evac=$OPTARG ;; g) eg=$OPTARG ;; O) out=$OPTARG ;;
               *) sed -n '2,30p' "$0"; exit 0 ;; esac
done
shift $((OPTIND-1))
[[ -f "$fc" ]] || die "give the ferrocene voltammogram with -F"
[[ "$dir" =~ ^(red|ox)$ ]] || die "-d must be red or ox"
mkdir -p "$out"
scale=$([[ $unit == mV ]] && echo 0.001 || echo 1)

# clean FILE -> "E(V) I(A) sweep" for the chosen cycle; sweep = +1 anodic, -1 cathodic
clean() {
    awk -v sc="$scale" -v want="$cyc" '
        { sub(/\r$/, ""); gsub(/,/, " ") }
        NF >= 2 && $1 ~ /^[-+]?[0-9.]+([eE][-+]?[0-9]+)?$/ && $2 ~ /^[-+]?[0-9.]+([eE][-+]?[0-9]+)?$/ {
            if ($1 + 0 == 0 && $2 + 0 == 0) next
            n++; E[n] = $1 * sc; I[n] = $2 }
        END {
            if (n < 10) { print "too few data points" > "/dev/stderr"; exit 1 }
            # direction of every step, smoothed over a few points to ignore noise
            for (i = 2; i <= n; i++) d[i] = (E[i] > E[i-1]) ? 1 : (E[i] < E[i-1] ? -1 : d[i-1]); d[1] = d[2]
            first = d[1]; c = 1; cyc[1] = 1; turns = 0
            for (i = 2; i <= n; i++) {
                if (d[i] != d[i-1]) { turns++; if (turns % 2 == 0 && d[i] == first) c++ }
                cyc[i] = c }
            nc = c; use = (want == "last") ? nc : want + 0
            # the last cycle is often incomplete: fall back to the one before
            cnt = 0; for (i = 1; i <= n; i++) if (cyc[i] == use) cnt++
            if (want == "last" && nc > 1 && cnt < 0.6 * n / nc) use = nc - 1
            for (i = 1; i <= n; i++) if (cyc[i] == use) printf "%.6f %.6e %d\n", E[i], I[i], d[i]
        }' "$1"
}

# ---------- ferrocene ----------
clean "$fc" > "$out/fc.dat" || die "cannot read $fc"
read -r epa epc e12 dep <<< "$(awk -v w="$fw" -v sc="$scale" '
    BEGIN { if (w != "") { split(w, a, ":"); lo = a[1] * sc; hi = a[2] * sc } }
    w == "" || ($1 >= lo && $1 <= hi) {
        if ($3 > 0 && (ia == "" || $2 > ia)) { ia = $2; ea = $1 }
        if ($3 < 0 && (ic == "" || $2 < ic)) { ic = $2; ec = $1 } }
    END { printf "%.4f %.4f %.4f %.1f\n", ea, ec, (ea + ec) / 2, (ea - ec) * 1000 }' "$out/fc.dat")"
info "Fc+/Fc: Epa = $epa V, Epc = $epc V, E1/2 = $e12 V, dEp = $dep mV"
awk -v d="$dep" 'BEGIN { if (d > 120 || d < 50) printf "warning: dEp = %s mV; a reversible one-electron couple gives about 59 mV at 25 C plus uncompensated resistance\n", d > "/dev/stderr" }'
echo "file,Epa_V,Epc_V,E12_V,dEp_mV" > "$out/ferrocene.csv"
echo "$(basename "$fc"),$epa,$epc,$e12,$dep" >> "$out/ferrocene.csv"

# ---------- analytes ----------
echo "file,wave,E_onset_vs_Fc_V,E_peak_vs_Fc_V,i_peak_A,E_LUMO_eV,E_HOMO_eV,Eg_opt_eV" > "$out/onsets.csv"
for f in "$@"; do
    n=$(strip_ext "$f")
    clean "$f" | awk -v s="$e12" '{ printf "%.6f %.6e %d\n", $1 - s, $2, $3 }' > "$out/$n-vsFc.dat" || { warn "$n: unreadable"; continue; }
    res=$(awk -v dir="$dir" -v ow="$ow" -v bw="$bw" '
        BEGIN { want = (dir == "red") ? -1 : 1; sgn = want
                if (ow != "") { split(ow, a, ":"); olo = a[1]; ohi = a[2] } else { olo = -1e9; ohi = 1e9 }
                if (bw != "") { split(bw, b, ":"); blo = b[1]; bhi = b[2] } }
        $3 == want && $1 >= olo && $1 <= ohi { n++; E[n] = $1; I[n] = $2 }
        function fit(i0, i1,   k, sx, sy, sxx, sxy, m) { sx = sy = sxx = sxy = 0; m = 0
            for (k = i0; k <= i1; k++) { sx += E[k]; sy += I[k]; sxx += E[k]^2; sxy += E[k]*I[k]; m++ }
            SL = (m * sxy - sx * sy) / (m * sxx - sx * sx); IC = (sy - SL * sx) / m }
        END {
            if (n < 20) { print "ERR too few points on the requested sweep"; exit }
            # peak: largest current in the direction of the wave (most negative for a reduction)
            ip = I[1]; ep = E[1]; jp = 1
            for (i = 1; i <= n; i++) if (sgn * I[i] > sgn * ip) { ip = I[i]; ep = E[i]; jp = i }
            # steepest rise before the peak, slope from a 5-point window
            best = 0; js = 0
            for (i = 3; i <= jp - 2; i++) { r = sgn * (I[i+2] - I[i-2]); if (r > best) { best = r; js = i } }
            if (js == 0) { print "ERR no rising current before the peak"; exit }
            fit(js - 2, js + 2); s1 = SL; c1 = IC
            # baseline: requested window or the first 15 % of the sweep
            if (bw != "") { b0 = 0; b1 = 0; for (i = 1; i <= n; i++) if (E[i] >= blo && E[i] <= bhi) { if (!b0) b0 = i; b1 = i } }
            else { b0 = 1; b1 = int(0.15 * n); if (b1 >= js - 2) b1 = js - 3 }
            if (b1 - b0 < 2) { print "ERR baseline window too small"; exit }
            fit(b0, b1); s0 = SL; c0 = IC
            eon = (c0 - c1) / (s1 - s0)
            printf "%.4f %.4f %.4e\n", eon, ep, ip }' "$out/$n-vsFc.dat")
    [[ "$res" == ERR* ]] && { warn "$n: ${res#ERR }"; continue; }
    read -r eon ep ip <<< "$res"
    awk -v f="$n" -v d="$dir" -v e="$eon" -v ep="$ep" -v ip="$ip" -v vac="$evac" -v g="$eg" 'BEGIN {
        if (d == "red") { lumo = -(e + vac); homo = (g != "") ? lumo - g : "" }
        else            { homo = -(e + vac); lumo = (g != "") ? homo + g : "" }
        printf "%s,%s,%.3f,%.3f,%s,%s,%s,%s\n", f, d, e, ep, ip, (lumo == "" ? "" : sprintf("%.2f", lumo)), (homo == "" ? "" : sprintf("%.2f", homo)), g }' >> "$out/onsets.csv"
    info "$n: $dir onset $eon V vs Fc+/Fc (peak $ep V)"
done

# ---------- plot ----------
if command -v gnuplot >/dev/null && (($# > 0)); then
    plots=(); for f in "$@"; do n=$(strip_ext "$f"); [[ -s "$out/$n-vsFc.dat" ]] && plots+=("'$out/$n-vsFc.dat' u 1:(\$2*1e6) w l lw 2 t '$n'"); done
    gnuplot <<EOG
set terminal pngcairo size 1000,700 enhanced font 'Arial,13'
set output '$out/cv-vs-fc.png'
set xlabel 'E (V vs Fc^+/Fc)'; set ylabel 'I ({/Symbol m}A)'; set grid lc rgb '#dddddd'; set key bottom right
plot $(IFS=,; echo "${plots[*]}")
EOG
    info "plot: $out/cv-vs-fc.png"
fi
