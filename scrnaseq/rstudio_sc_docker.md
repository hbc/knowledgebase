# Docker image with rstudio for single cell analysis

## Description

Docker images for single cell analysis.

All docker images contain an rstudio installation with some helpful packages for singlecell analysis. It also includes a conda environment to deal with necessary python packages (like umap-learn).

Docker Rstudio images are obtained from [rocker/rstudio](https://hub.docker.com/r/rocker/rstudio).

## R version and Bioconductor

The R and Bioconductor versions are specified in the image name (along with the OS version):

Example:
`singlecell-base:R.4.0.3-BioC.3.11-ubuntu_20.04`

## Use 

`docker run -d -p 8787:8787 --name <container_name> -e USER='rstudio' -e PASSWORD='rstudioSC' -e ROOT=TRUE -v <host_folder>:/home/rstudio/projects vbarrerab/<docker_image>)`

`-e DISABLE_AUTH=true` option can be added to avoid Rstudio login prompt. Only use on local machine.

This instruction will download and launch a container using the singlecell image. Once launch, it can be access through a web browser with the URL 8787:8787 or localhost:8787.

### Important parameters

* -v option is mounting a folder from the host in the container. This allows for data transfer between the host and the container. **This can only be done when creating the container!**

* --name assigns a name to the container. Helpful to keep thins tidy.
* -e ROOT=TRUE options provides root access, in case more tweaking inside the container is necessary.
* -p 8787:<port> Change the local port to access the container. **This can only be done when creating the container!**
* FYI: The working directory will be set as /home/rstudio, not /home/rstudio/projects as default behavior.

## Resources

The dockerfile and other configuration files can be found on:

https://github.com/vbarrera/docker_configuration

The docker images: 

vbarrerab/singlecell-base

## Available images:

- R.4.0.2-BioC.3.11-ubuntu_20.04
- R.4.0.3-BioC.3.11-ubuntu_20.04

**Important:**

Docker changed its policies to only keep images that have been modified in the last 6 months. This means that previous images will eventually disappear. For previous versions. Check with availability with @vbarrera.

# Bibliography

Inspired by:

https://www.r-bloggers.com/running-your-r-script-in-docker/
  
# Other resources
Using Singularity Containers on the Cluster: https://docs.rc.fas.harvard.edu/kb/singularity-on-the-cluster/
