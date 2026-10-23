#!/usr/bin/env bash
# =============================================================================
#  uv-onset.sh - absorption onset and optical gap from UV-Vis spectra
#
#  usage: uv-onset.sh [-w LO:HI] [-f 0.1] [-O outdir] spectrum1.txt ...
#    spectra: two columns, wavelength (nm) and absorbance (any header skipped)
#    -w   wavelength window to search (default: whole spectrum)
#    -f   a band counts as the lowest energy band when its maximum is at least
#         this fraction of the strongest absorbance (default 0.1)
#  lambda_max  = maximum of the longest wavelength band
#  lambda_onset = crossing of the tangent at the steepest point of its red edge
#                 with the baseline (fitted to the last 10 % of the window)
#  Eg_opt = 1239.84 / lambda_onset (eV)
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"
win=""; frac=0.1; out=uv-results
while getopts "w:f:O:h" o; do case $o in w) win=$OPTARG ;; f) frac=$OPTARG ;; O) out=$OPTARG ;; *) sed -n '2,15p' "$0"; exit 0 ;; esac; done
shift $((OPTIND-1)); (($#)) || die "give spectra"
mkdir -p "$out"; csv="$out/uv-onsets.csv"; echo "file,lambda_max_nm,A_max,lambda_onset_nm,Eg_opt_eV" > "$csv"
for f in "$@"; do
    awk -v win="$win" -v fr="$frac" -v name="$(strip_ext "$f")" '
        BEGIN { if (win != "") { split(win, w, ":"); lo = w[1]; hi = w[2] } else { lo = -1e9; hi = 1e9 } }
        { sub(/\r$/, ""); gsub(/,/, " ") }
        $1 ~ /^[0-9.]+$/ && $2 ~ /^[-+]?[0-9.]+([eE][-+]?[0-9]+)?$/ && $1 >= lo && $1 <= hi { n++; L[n] = $1; A[n] = $2 }
        END {
            if (n < 10) { print "too few points in " name > "/dev/stderr"; exit 1 }
            # sort by wavelength (spectrometers often export from long to short)
            for (i = 2; i <= n; i++) { l = L[i]; a = A[i]; j = i - 1; while (j > 0 && L[j] > l) { L[j+1] = L[j]; A[j+1] = A[j]; j-- } L[j+1] = l; A[j+1] = a }
            gmax = 0; for (i = 1; i <= n; i++) if (A[i] > gmax) gmax = A[i]
            # longest wavelength local maximum above the threshold (3-point check)
            jm = 0; for (i = n - 1; i >= 2; i--) if (A[i] >= A[i-1] && A[i] >= A[i+1] && A[i] >= fr * gmax) { jm = i; break }
            if (!jm) { print "no band found in " name > "/dev/stderr"; exit 1 }
            # steepest descent on the red side
            best = 0; js = 0; for (i = jm + 2; i <= n - 2; i++) { r = (A[i-2] - A[i+2]) / (L[i+2] - L[i-2]); if (r > best) { best = r; js = i } }
            if (!js) { print "no red edge in " name > "/dev/stderr"; exit 1 }
            sx = sy = sxx = sxy = 0; for (k = js - 2; k <= js + 2; k++) { sx += L[k]; sy += A[k]; sxx += L[k]^2; sxy += L[k]*A[k] }
            s1 = (5 * sxy - sx * sy) / (5 * sxx - sx * sx); c1 = (sy - s1 * sx) / 5
            b0 = int(0.9 * n); if (b0 <= js + 2) b0 = js + 3
            m = sx = sy = sxx = sxy = 0; for (k = b0; k <= n; k++) { m++; sx += L[k]; sy += A[k]; sxx += L[k]^2; sxy += L[k]*A[k] }
            if (m >= 3) { s0 = (m * sxy - sx * sy) / (m * sxx - sx * sx); c0 = (sy - s0 * sx) / m } else { s0 = 0; c0 = 0 }
            lon = (c0 - c1) / (s1 - s0)
            printf "%s,%.1f,%.4f,%.1f,%.3f\n", name, L[jm], A[jm], lon, 1239.84 / lon }' "$f" >> "$csv" || warn "$f: no onset"
done
info "UV onsets: $csv"; cat "$csv" >&2
