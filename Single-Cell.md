# single cell RNAseq

## Shared installation in O2

There is a Rlib shared installation in O2 with seurat and zinbvawe and tidyverse packages.

To use it, load:

```
module load module load gcc/6.2.0 R/3.4.1 gsl
```

Add this to `~/.Renviron`

```
R_LIBS_USER=/n/data1/cores/bcbio/R/library/3.4-bioc-release/library
R_MAX_NUM_DLLS=150
```

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