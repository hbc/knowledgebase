---
title: Single Cell Installations
description: This code helps with installing tools and setting up docker for single cell rnaseq.
category: research
subcategory: scrnaseq
tags: [hpc, R]
---


# single cell RNAseq

## Shared installation in O2

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

## Docker container with bcbioSingleCell

There is a docker container with the last version of R and Seurat installed and kind of updated bcbioSingleCell package.

To set up the Docker container to run the single-cell RNA-Seq analysis using bcbioSingleCell and Seurat, we need to perform a few steps. 

> **NOTE:** There are nice tutorials online to help understand the basics of Docker, including [Docker for beginners](https://docker-curriculum.com) and [Titus' materials](http://angus.readthedocs.io/en/2016/week3/CTB_docker.html).

To set up Docker on a local computer, install the Docker application, then run the `docker pull` command using the `r3.5-bsc0.1.5` tag available from `lpantano/bcbiosinglecell`(available tags can be found [here](https://hub.docker.com/r/lpantano/bcbiosinglecell/tags/)):

```bash
docker pull lpantano/bcbiosinglecell:r3.5-bsc0.1.5
```

After pulling the docker image, run the container as detailed at https://hub.docker.com/r/lpantano/bcbiosinglecell/. This site also has explanations on how to install/update packages when you have it running in your local computer.

# zinbwave and DESeq2
If you have a complicated experimental design, zinbwave can assign invalid weights; DESeq2 will fail with a `weights.ok not all TRUE` message. Mike and I fixed that in the devel version of `DESeq2` which you can install with `devtools::install_github("mikelove/DESeq2")`.
