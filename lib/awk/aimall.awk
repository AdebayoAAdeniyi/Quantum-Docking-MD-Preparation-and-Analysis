# =============================================================================
#  aimall.awk - bond critical points from AIMAll .sum or .mgpviz files
#  usage: awk [-v name=NAME] [-v header=1] -f aimall.awk file.sum
#  prints name,cp,atom_A,atom_B,rho,DelSqRho,ellipticity,G,V,H,minusG_over_V,bond_class
#   H = G + V (total energy density, Cremer and Kraka 1984)
#   bond_class from the signs of DelSqRho and H and the value of -G/V
#   (Espinosa et al., J. Chem. Phys. 2002, 117, 5529):
#     shared   DelSqRho < 0                       (covalent)
#     interm.  DelSqRho > 0, H < 0                (partly covalent: dative, metal-ligand)
#     closed   DelSqRho > 0, H > 0                (ionic, hydrogen bond, van der Waals)
#  Values in atomic units, as AIMAll writes them.
# =============================================================================
function num(s) { gsub(/[dD]/, "E", s); return s + 0 }
function flush() {
    if (!inb) return
    H = G + V; r = (V != 0 ? -G / V : 0)
    cls = (L < 0 ? "shared" : (H < 0 ? "intermediate" : "closed-shell"))
    printf "%s,%s,%s,%s,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.4f,%s\n", name, cp, a, b, rho, L, ell, G, V, H, r, cls
    inb = 0
}
BEGIN { if (name == "") name = "mol"; if (header) print "name,cp,atom_A,atom_B,rho,DelSqRho,ellipticity,G,V,H,minusG_over_V,bond_class" }
/^CP# +[0-9]+/ { flush(); cp = $2; next }
/Type = \(3,-1\) BCP/ { inb = 1; a = $(NF - 1); b = $NF; rho = L = ell = G = V = ""; next }
/Type = \(3,[-+]?[13]\)/ && !/BCP/ { flush(); next }
inb && /^ *Rho = /              { rho = num($3) }
inb && /^ *DelSqRho = /         { L = num($3) }
inb && /^ *Bond Ellipticity = / { ell = num($4) }
inb && /^ *G = /                { G = num($3) }
inb && /^ *V = /                { V = num($3) }
END { flush() }
