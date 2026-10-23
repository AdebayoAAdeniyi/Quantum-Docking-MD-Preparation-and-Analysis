#!/usr/bin/env bash
# =============================================================================
#  build-bead-spring.sh - coarse-grained polymer melt for testing the LAMMPS
#                         property workflows (reduced Lennard-Jones units)
#
#  Chains of N beads (random walks with bond length 0.97 sigma) at number density
#  rho, harmonic bonds, LJ 12-6 with a 2.5 sigma cutoff (attractive, so the melt
#  vitrifies on cooling). Overlaps from the random start are removed in the
#  equilibration by a soft-core push-off before LJ is switched on.
#
#  usage: build-bead-spring.sh [-c 50] [-n 20] [-r 0.85] [-o polymer]
#  writes polymer/system.data and polymer/system.in (styles and coefficients)
# =============================================================================
source "$(dirname "$0")/../../lib/common.sh"
nc=50; nb=20; rho=0.85; out=polymer; seed=12345
while getopts "c:n:r:o:s:h" o; do case $o in c) nc=$OPTARG ;; n) nb=$OPTARG ;; r) rho=$OPTARG ;; o) out=$OPTARG ;; s) seed=$OPTARG ;; *) sed -n '2,13p' "$0"; exit 0 ;; esac; done
mkdir -p "$out"
awk -v nc="$nc" -v nb="$nb" -v rho="$rho" -v seed="$seed" 'BEGIN {
    srand(seed); N = nc * nb; L = (N / rho) ^ (1/3); pi = 3.141592653589793
    printf "bead-spring melt: %d chains x %d beads, rho = %s\n\n%d atoms\n%d bonds\n\n1 atom types\n1 bond types\n\n", nc, nb, rho, N, nc * (nb - 1)
    printf "0 %.6f xlo xhi\n0 %.6f ylo yhi\n0 %.6f zlo zhi\n\nMasses\n\n1 1.0\n\nAtoms # bond\n\n", L, L, L
    a = 0
    for (c = 1; c <= nc; c++) {
        x = rand() * L; y = rand() * L; z = rand() * L
        for (b = 1; b <= nb; b++) {
            if (b > 1) { do { ct = 2 * rand() - 1; ph = 2 * pi * rand(); st = sqrt(1 - ct * ct)
                              dx = 0.97 * st * cos(ph); dy = 0.97 * st * sin(ph); dz = 0.97 * ct
                              # avoid immediate back-folding onto the bead before last
                          } while (b > 2 && (x + dx - px) ^ 2 + (y + dy - py) ^ 2 + (z + dz - pz) ^ 2 < 1.0)
                          px = x; py = y; pz = z; x += dx; y += dy; z += dz }
            a++; printf "%d %d 1 %.5f %.5f %.5f\n", a, c, x, y, z } }
    printf "\nBonds\n\n"; k = 0
    for (c = 0; c < nc; c++) for (b = 1; b < nb; b++) { k++; printf "%d 1 %d %d\n", k, c * nb + b, c * nb + b + 1 } }' > "$out/system.data"
cat > "$out/system.in" <<'EOT'
# bead-spring polymer, reduced LJ units
units           lj
atom_style      bond
boundary        p p p
read_data       ${data}
bond_style      harmonic
bond_coeff      1 500.0 0.97
special_bonds   lj 0.0 1.0 1.0
pair_style      lj/cut 2.5
pair_coeff      1 1 1.0 1.0 2.5
pair_modify     shift yes
neighbor        0.4 bin
neigh_modify    every 1 delay 0 check yes
EOT
info "$out/system.data: $((nc * nb)) beads in $nc chains; $out/system.in"
