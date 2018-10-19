---
title: Single Cell Installations
description: This code helps with installing tools and setting up docker for single cell rnaseq.
category: research
subcategory: scrnaseq
tags: [hpc, R]
---


# single cell RNAseq

## R Shared installation in O2

There is a Rlib shared installation in O2 with seurat and zinbvawe and tidyverse packages.

To use it, load:

```
module load gcc/6.2.0 R/3.4.1 gsl hdf5
```

Add this to `~/.Renviron`

```
R_LIBS_USER=/n/data1/cores/bcbio/R/library/3.4-bioc-release/library
R_MAX_NUM_DLLS=150
```

There is also a Rlib installation in O2 with seurat and zinbvawe and tidyverse packages for R.3.5.1


To use it, load:

```
module load gcc/6.2.0 R/3.5.1 gsl hdf5
```

Add this to `~/.Renviron`

```
R_LIBS_USER=/n/data1/cores/bcbio/R/library/3.5-bioc-release/library
R_MAX_NUM_DLLS=150
```

### Notes:

During Zinbwave installation an error may arise with the gsl package. To solve it:

With gsl module loaded find which flags are set with: 
```
gsl-config --cflags

gsl-config --libs
```
Export them before loading R as (the indicated flags are O2's at 08/02/2018):

```
export CFLAGS="-I/n/app/gsl/2.3/include"; export LDFLAGS="-L/n/app/gsl/2.3/lib -lgsl -lgslcblas -lm"
```

Now you can install the package normally.

(source: https://stackoverflow.com/questions/24781125/installing-r-gsl-package-on-mac)

Then just type `R` or use `Rscript`.

## How to run main analysis

Read the full documentation at: https://github.com/marypiper/bcbio_single_cell_workflow/tree/master/lessons. It provides with step by step explanation on how to analyze single cell data inside HBC core. If you see something out of date, please open a PR at that repository.

## Visualization

Reports to visualize markers and clusters can be found here: http://bioinformatics.sph.harvard.edu/hbcABC/articles/singlecell-tutorial.html. It includes simple QC report and clustering and markers reports.

# zinbwave and DESeq2
If you have a complicated experimental design, zinbwave can assign invalid weights; DESeq2 will fail with a `weights.ok not all TRUE` message. Mike and I fixed that in the devel version of `DESeq2` which you can install with `devtools::install_github("mikelove/DESeq2")`.

# Use DESeq2 only for UMI disambiguated protocols
zinbwave will cause you to miss many of the low expressors that would be good markers. See 
https://github.com/roryk/zinbwave-deseq2-indrop and https://support.bioconductor.org/p/112163/
