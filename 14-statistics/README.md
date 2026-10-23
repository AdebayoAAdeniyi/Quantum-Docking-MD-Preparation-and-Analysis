# 14 · Statistics and trajectory PCA (R)

Two R scripts that close the loop: one relates computed properties to each other and to measured quantities, the other analyses the collective motion in an MD trajectory.

```
14-statistics/
├── r-analysis/property-stats.R   correlations, PCA (FactoMineR), multiple linear regression
└── bio3d/bio3d-analysis.R        trajectory PCA, RMSF, dynamic cross-correlation (bio3d)
```

## Property statistics

```bash
Rscript 14-statistics/r-analysis/property-stats.R -i properties.csv -y logk_exp \
        -x "gap_eV,omega_eV,dipole_D" -o stats-results
```

Takes any CSV with one row per molecule, such as the tables the extract scripts produce. It writes the correlation matrix (CSV and a labelled heat map), a PCA with scores, loadings, explained variance and a biplot, and, when a response column is given with `-y`, the single-variable correlations ranked by |r| and a multiple linear regression with its coefficients, R², adjusted R² and an observed-versus-predicted plot.

It warns when there are as many predictors as molecules, because R² is then meaningless. That is an easy trap when a series has five or six compounds, which is common in this kind of work.

Run on the five benzophenones screened with xTB in module 04, PC1 and PC2 carry 58 % and 31 % of the variance of ten computed properties.

## Trajectory analysis

```bash
Rscript 14-statistics/bio3d/bio3d-analysis.R -p solute.pdb -t prod.nc -s calpha -k 3
```

Fits the frames, then writes per-residue RMSF, the principal components of the Cartesian coordinates with their eigenvalues and the projection of every frame, the per-residue contribution to PC1 and PC2, a clustering of the frames in the PC1/PC2 plane, the dynamic cross-correlation matrix, and `pc1.pdb`, a trajectory along PC1 for viewing the dominant motion.

## Corrections made to the original scripts

* The old R scripts had the file names, column numbers and plot ranges typed into them, so each new data set meant editing the script. Everything is now a command line option and columns are chosen by name.
* The PCA, correlation and regression steps were separate scripts that each re-read and re-cleaned the data; they are one pass now, and constant columns (which make the correlation matrix undefined) are dropped with a message.
* bio3d, DSSP and the R tarballs that were kept beside the original scripts are not redistributed; they install from CRAN and conda-forge (see `install/README.md`).

## Tested

`tests/test-14-statistics.sh` runs both scripts: the property statistics on the real xTB screening table of module 04 (with a synthetic response column, so the regression path is exercised), and bio3d on a short trajectory. Both need R with `FactoMineR` and `bio3d`.
