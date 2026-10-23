#!/usr/bin/env Rscript
# =============================================================================
#  bio3d-analysis.R - principal component analysis and dynamic cross-correlation
#                     of an MD trajectory (bio3d)
#
#  usage: Rscript bio3d-analysis.R -p topology.pdb -t trajectory.{dcd,nc,xtc}
#                                  [-s "calpha"] [-k 3] [-o bio3d-results]
#    -p  a PDB of the solute (cpptraj: trajout x.pdb onlyframes 1; gmx trjconv -dump 0)
#    -t  trajectory in a format bio3d reads (DCD or NetCDF; convert XTC first)
#    -s  atom selection (default "calpha")
#    -k  number of conformational clusters along PC1/PC2 (default 3)
#  Outputs: rmsf.csv, pca-eigenvalues.csv, pca-projection.csv (frames on PC1-PC3),
#  pc1-loadings.csv (contribution of each residue), dccm.csv (cross-correlation),
#  pca-clusters.csv, and PNG plots; pc1.pdb is the PC1 trajectory for viewing.
#  Replaces run-bio3d / bio3d-full.input and the separate network script.
# =============================================================================
suppressMessages(library(bio3d))
args <- commandArgs(trailingOnly = TRUE)
ga <- function(f, d = NA) { i <- match(f, args); if (is.na(i) || i == length(args)) d else args[i + 1] }
pdbf <- ga("-p"); trjf <- ga("-t"); sel <- ga("-s", "calpha"); k <- as.integer(ga("-k", 3)); out <- ga("-o", "bio3d-results")
if (is.na(pdbf) || is.na(trjf)) stop("give -p topology.pdb and -t trajectory")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
pdb <- read.pdb(pdbf)
trj <- if (grepl("\\.nc$|\\.ncdf$", trjf)) read.ncdf(trjf) else read.dcd(trjf)
ca <- atom.select(pdb, sel)
cat(sprintf("%d frames, %d selected atoms\n", nrow(trj), length(ca$xyz) / 3))
xyz <- fit.xyz(fixed = pdb$xyz, mobile = trj, fixed.inds = ca$xyz, mobile.inds = ca$xyz)
resno <- pdb$atom$resno[ca$atom]
# RMSF
r <- rmsf(xyz[, ca$xyz])
write.csv(data.frame(residue = resno, rmsf = round(r, 4)), file.path(out, "rmsf.csv"), row.names = FALSE)
# PCA of the Cartesian coordinates
pc <- pca.xyz(xyz[, ca$xyz])
write.csv(data.frame(PC = seq_along(pc$L), eigenvalue = round(pc$L, 5),
                     variance_percent = round(100 * pc$L / sum(pc$L), 2)), file.path(out, "pca-eigenvalues.csv"), row.names = FALSE)
write.csv(round(pc$z[, 1:3], 4), file.path(out, "pca-projection.csv"))
write.csv(data.frame(residue = resno, PC1 = round(pc$au[, 1], 4), PC2 = round(pc$au[, 2], 4)),
          file.path(out, "pc1-loadings.csv"), row.names = FALSE)
cat(sprintf("PC1 %.1f%%, PC2 %.1f%% of the fluctuation\n", 100 * pc$L[1] / sum(pc$L), 100 * pc$L[2] / sum(pc$L)))
# clusters in the PC1-PC2 plane
grp <- cutree(hclust(dist(pc$z[, 1:2])), k = k)
write.csv(data.frame(frame = seq_along(grp), cluster = grp), file.path(out, "pca-clusters.csv"), row.names = FALSE)
# dynamic cross-correlation
cij <- dccm(xyz[, ca$xyz])
write.csv(round(cij, 3), file.path(out, "dccm.csv"))
png(file.path(out, "pca-projection.png"), width = 900, height = 450)
par(mfrow = c(1, 2))
plot(pc$z[, 1], pc$z[, 2], col = grp, pch = 19, xlab = "PC1 (A)", ylab = "PC2 (A)", main = "conformational space")
plot(100 * pc$L / sum(pc$L), type = "h", lwd = 6, xlim = c(0, 10), xlab = "PC", ylab = "% variance", main = "variance per PC")
invisible(dev.off())
png(file.path(out, "rmsf.png"), width = 900, height = 450)
plot(resno, r, type = "l", lwd = 2, col = "#1f5f99", xlab = "residue", ylab = "RMSF (A)")
invisible(dev.off())
png(file.path(out, "dccm.png"), width = 800, height = 750)
plot(cij, resno = pdb, main = "dynamic cross-correlation")
invisible(dev.off())
mktrj(pc, pc = 1, file = file.path(out, "pc1.pdb"))
cat(sprintf("results in %s\n", out))
