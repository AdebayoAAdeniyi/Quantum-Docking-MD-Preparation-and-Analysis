#!/usr/bin/env Rscript
# =============================================================================
#  property-stats.R - correlation, PCA and linear models for a table of computed
#                     (or measured) molecular properties
#
#  usage: Rscript property-stats.R -i properties.csv [-n name] [-y response]
#                                  [-x "col1,col2"] [-o outdir] [-s]
#    -i  CSV with one row per molecule
#    -n  column holding the molecule name (default: the first column)
#    -y  response variable for the linear model and the correlation ranking
#        (e.g. an experimental activity or rate); omit for PCA and correlations only
#    -x  predictors (default: every other numeric column)
#    -s  scale the variables to unit variance before PCA (recommended when the
#        properties have different units; default TRUE unless -s is given twice)
#  Outputs: correlation-matrix.csv and .png, pca-scores.csv, pca-loadings.csv,
#  pca-variance.csv, pca-biplot.png, and with -y: regression.txt (coefficients,
#  R2, adjusted R2, F test), regression-vs-observed.png, correlations-with-y.csv.
#  Replaces the ad hoc R scripts (run-mult-Var-anal, plotwithLines.R, PLS/plspm).
# =============================================================================
suppressMessages({ library(FactoMineR) })
args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NA) { i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1] }
infile <- getarg("-i"); if (is.na(infile)) stop("give -i properties.csv")
outdir <- getarg("-o", "stats-results"); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dat <- read.csv(infile, stringsAsFactors = FALSE, check.names = TRUE)
namecol <- getarg("-n", names(dat)[1])
rownames(dat) <- make.unique(as.character(dat[[namecol]]))
y <- getarg("-y", NA)
num <- dat[, sapply(dat, is.numeric), drop = FALSE]
xs <- getarg("-x", NA)
if (!is.na(xs)) num <- num[, intersect(unlist(strsplit(xs, ",")), names(num)), drop = FALSE]
num <- num[, apply(num, 2, function(v) sd(v, na.rm = TRUE) > 0), drop = FALSE]
if (ncol(num) < 2) stop("fewer than two numeric columns with any variation")
cat(sprintf("%d molecules, %d numeric variables\n", nrow(num), ncol(num)))
# ---- correlations ----
cm <- cor(num, use = "pairwise.complete.obs")
write.csv(round(cm, 4), file.path(outdir, "correlation-matrix.csv"))
png(file.path(outdir, "correlation-matrix.png"), width = 900, height = 800)
par(mar = c(9, 9, 2, 2))
image(1:ncol(cm), 1:ncol(cm), t(cm[ncol(cm):1, ]), col = colorRampPalette(c("#2166ac", "white", "#b2182b"))(41),
      zlim = c(-1, 1), axes = FALSE, xlab = "", ylab = "")
axis(1, 1:ncol(cm), colnames(cm), las = 2, cex.axis = 0.8); axis(2, 1:ncol(cm), rev(colnames(cm)), las = 2, cex.axis = 0.8)
for (i in 1:ncol(cm)) for (j in 1:ncol(cm)) text(i, ncol(cm) - j + 1, sprintf("%.2f", cm[j, i]), cex = 0.6)
invisible(dev.off())
# ---- PCA ----
scale.it <- !(sum(args == "-s") >= 2)
p <- PCA(num, scale.unit = scale.it, graph = FALSE)
write.csv(round(p$ind$coord, 4), file.path(outdir, "pca-scores.csv"))
write.csv(round(p$var$coord, 4), file.path(outdir, "pca-loadings.csv"))
write.csv(round(p$eig, 4), file.path(outdir, "pca-variance.csv"))
png(file.path(outdir, "pca-biplot.png"), width = 1000, height = 500)
par(mfrow = c(1, 2))
plot(p, choix = "ind", title = "PCA scores")
plot(p, choix = "var", title = "PCA loadings")
invisible(dev.off())
cat(sprintf("PC1 %.1f%%, PC2 %.1f%% of the variance\n", p$eig[1, 2], p$eig[2, 2]))
# ---- regression ----
if (!is.na(y) && y %in% names(dat)) {
    preds <- setdiff(names(num), y)
    r <- sapply(preds, function(v) cor(num[[v]], dat[[y]], use = "pairwise.complete.obs"))
    r <- sort(r[!is.na(r)], decreasing = TRUE)
    write.csv(data.frame(variable = names(r), r = round(r, 4), r2 = round(r^2, 4)),
              file.path(outdir, "correlations-with-y.csv"), row.names = FALSE)
    f <- as.formula(paste(y, "~", paste(preds, collapse = " + ")))
    m <- lm(f, data = cbind(num, setNames(data.frame(dat[[y]]), y)))
    capture.output(summary(m), file = file.path(outdir, "regression.txt"))
    s <- summary(m)
    cat(sprintf("linear model on %d predictors: R2 = %.4f, adjusted R2 = %.4f\n", length(preds), s$r.squared, s$adj.r.squared))
    if (nrow(num) <= length(preds) + 1) cat("warning: as many predictors as molecules; R2 is meaningless here, use fewer variables\n")
    png(file.path(outdir, "regression-vs-observed.png"), width = 700, height = 650)
    obs <- dat[[y]]; prd <- predict(m)
    plot(obs, prd, pch = 19, col = "#1f5f99", xlab = paste("observed", y), ylab = paste("predicted", y),
         main = sprintf("R2 = %.3f", s$r.squared)); abline(0, 1, lty = 2)
    text(obs, prd, rownames(num), pos = 3, cex = 0.7)
    invisible(dev.off())
}
cat(sprintf("results in %s\n", outdir))
