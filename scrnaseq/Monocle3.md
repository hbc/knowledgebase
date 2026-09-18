# Installation tips for Monocle3 on O2

## You need to have certain modules loaded in your terminal or O2 Portal session

1. Load required modules

- `udunits` is required to install the `units` R package
- `gdal`, `geos`, and, `proj` are required to install the `sf` R package

Terminal:
```
ml gcc/14.2.0
ml udunits/2.2.28 gdal/3.11.3 geos/3.13.1 proj/9.6.2 R/4.4.2
```

Package list to request O2 Portal RStudio session:
```
gcc/14.2.0 udunits/2.2.28 gdal/3.11.3 geos/3.13.1 proj/9.6.2 R/4.4.2
```

2. Start R

3. Define your personal R package library

We recommend making a separate library for each version of R.

```
userlib <- "~/R/4.4.2"
```

4. Install `units`

```
install.packages("units", lib = userlib, type = "source", repos = "https://cloud.r-project.org")
```

5. Install `sf`

```
install.packages("sf",    lib = userlib, type = "source", repos = "https://cloud.r-project.org")
```

6. Install `spdep`

```
install.packages("spdep", lib = userlib, type = "source", repos = "https://cloud.r-project.org")

```

7. Install the packages that [Monocle3](https://cole-trapnell-lab.github.io/monocle3/docs/installation/) actually tells you that it needs
```
remotes::install_github('satijalab/seurat-wrappers')
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats', 'limma',
                       'lme4', 'S4Vectors', 'SingleCellExperiment', 'SummarizedExperiment',
                       'batchelor', 'HDF5Array', 'ggrastr'),
                     lib = userlib, lib.loc = userlib)
remotes::install_github("bnprks/BPCells/r")
remotes::install_github("cole-trapnell-lab/monocle3", lib = userlib, dependencies = FALSE)
```
