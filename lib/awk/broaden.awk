# =============================================================================
#  broaden.awk - turn stick spectra into smooth curves
#
#  input : lines "position intensity" (whitespace or comma separated)
#  usage : awk -v shape=lorentz|gauss -v width=W -v xmin=A -v xmax=B
#              [-v npts=2000] [-v scale=1] [-v uv=1] -f broaden.awk sticks
#    width   full width at half maximum, same unit as the positions
#    scale   factor applied to positions (frequency scaling, e.g. 0.967)
#    uv=1    positions in nm, intensities are oscillator strengths: the curve
#            is the molar absorptivity (L mol-1 cm-1) computed on an energy
#            grid with Gaussian bands, as in the GaussView/Gaussian white paper
#            eps(v) = 1.3062974e8 * f / sigma * exp(-((v - v_i) / sigma)^2),
#            v and sigma in cm-1; here width is sigma in eV (default 0.4 eV)
#  output: "x y" lines (x in nm for uv=1)
# =============================================================================
BEGIN { FS = "[ \t,]+"; if (npts == "") npts = 2000; if (scale == "") scale = 1; n = 0 }
/^[ \t]*#/ || NF < 2 { next }
{ i = ($1 == "" ? 2 : 1); if ($(i) + 0 == 0 && $(i) !~ /^-?0/) next
  n++; x[n] = $(i) * scale; y[n] = $(i + 1) }
END {
    pi = 3.141592653589793
    if (uv) {
        if (width == "") width = 0.4
        sig = width * 8065.544                 # eV -> cm-1
        for (k = 0; k < npts; k++) {
            lam = xmin + (xmax - xmin) * k / (npts - 1); v = 1e7 / lam; e = 0
            for (j = 1; j <= n; j++) { d = (v - 1e7 / x[j]) / sig; if (d * d < 50) e += 1.3062974e8 * y[j] / sig * exp(-d * d) }
            printf "%.3f %.4f\n", lam, e
        }
        exit
    }
    hw = width / 2; c = sqrt(4 * log(2) / pi) / width; a = 4 * log(2) / (width * width)
    for (k = 0; k < npts; k++) {
        xx = xmin + (xmax - xmin) * k / (npts - 1); s = 0
        for (j = 1; j <= n; j++) {
            d = xx - x[j]
            if (shape == "gauss") s += y[j] * c * exp(-a * d * d)
            else                  s += y[j] * (hw / pi) / (d * d + hw * hw)
        }
        printf "%.3f %.6g\n", xx, s
    }
}
