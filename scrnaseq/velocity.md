Velocity analysis is trajectory analysis based on spliced/unspliced RNA ratio.

It is quite popular https://www.nature.com/articles/s41586-018-0414-6,  
however, the original pipeline is not well supported:
https://github.com/velocyto-team/velocyto.R/issues

There is a new one from kallisto team:
https://bustools.github.io/BUS_notebooks_R/velocity.html

- module load gcc/6.2.0
- installed R-devel: https://www.r-bloggers.com/r-devel-in-parallel-to-regular-r-installation/  
because one of the packages wanted R4.0
- configure R with `./configure --enable-R-shlib` for rstudio
- remove conda from PATH to avoid using its libcurl
- module load boost/1.62.0
- module load hdf5/1.10.1
- installing velocyto.R: https://github.com/velocyto-team/velocyto.R/issues/86
