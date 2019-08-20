---
title: Docker image with rstudio for single cell analysis
description: Docker image with rstudio prepared for SC analysis
category: research
subcategory: scrnaseq
tag: [docker,single_cell,resource,R]
---

# Docker image with rstudio for single cell analysis

## Description

This docker image contains an rstudio installation with some helpful packages for singlecell analysis. It also includes a conda environment to deal with necessary pythong packages (like umap-learn).

The image contains the latest R version that comes with the latest rocker/rstudio image. In order to update it, it's necessary to rebuild the image. 

## Use 

`docker run -d -p 8787:8787 --name <container_name> -e USER='rstudio' -e PASSWORD='rstudioSC' -e ROOT=TRUE -v <host_folder>:/home/rstudio/projects vbarrerab/rstudio_singlecell:latest`

This instruction will download and launch a container using the rstudio_singlecell image. Once launch, it can be access through a webbrowser with the URL 8787:8787 or localhost:8787

### Important parameters

* -v option is mounting a folder from the host in the container. This allows for data transfer between the host and the container. **This can only be done when creating the container!**

* --name assigns a name to the container. Helpful to keep thins tidy.
* -e ROOT=TRUE options provides root access, in case more tweaking inside the container is necessary.

## Installed R packages:

* tidyverse
* tibble
* RCurl
* cowplot 
* SingleCellExperiment
* scater
* reticulate
* AnnotationHub
* ensembldb
* Seurat 
* rio
* devtools
* dplyr

## Resources

The dockerfile and other configuration files can be found on:

https://code.harvard.edu/vib406/docker_singlecell

The docker image: 

vbarrerab/rstudio_singlecell:latest

Inspired by:

https://www.r-bloggers.com/running-your-r-script-in-docker/
