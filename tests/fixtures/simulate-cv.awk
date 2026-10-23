# Writes a simulated cyclic voltammogram in the layout of the potentiostat text
# export used in the original work (header block, "Segment" lines, potential in
# mV, current in A). Shapes are idealised and only test the processing:
#   type=fc       reversible couple, E1/2 = e12 mV, peaks at E1/2 +/- 35 mV
#   type=red      reduction wave i = -ilim / (1 + exp((E - e0)/s)); its tangent
#                 onset is exactly e0 + 2 s
# usage: awk -v type=fc -v e12=450 -f simulate-cv.awk  (or type=red e0=-1700 s=30)
BEGIN {
    srand(7); for (h = 1; h <= 26; h++) print "Header line " h " (instrument settings)"
    cap = 2e-7; step = 2
    if (type == "fc") { a = -200; b = 900 } else { a = -2300; b = 300 }
    for (cyc = 1; cyc <= 2; cyc++) {
        print "Segment " (2 * cyc - 1)
        for (E = (type == "fc" ? a : b); (type == "fc" ? E <= b : E >= a); E += (type == "fc" ? step : -step)) out(E, type == "fc" ? 1 : -1)
        print "Segment " (2 * cyc)
        for (E = (type == "fc" ? b : a); (type == "fc" ? E >= a : E <= b); E += (type == "fc" ? -step : step)) out(E, type == "fc" ? -1 : 1)
    }
}
function out(E, dir,   i) {
    if (type == "fc") i = dir > 0 ? 1e-5 * exp(-((E - e12 - 35) / 40) ^ 2) + cap : -1e-5 * exp(-((E - e12 + 35) / 40) ^ 2) - cap
    else i = dir < 0 ? -2e-5 / (1 + exp((E - e0) / s)) - cap : -1e-5 / (1 + exp((E - e0) / s)) + cap
    printf "%.1f\t%.6E\n", E, i
}
