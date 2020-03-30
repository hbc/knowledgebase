RNA Velocity analysis is a trajectory analysis based on spliced/unspliced RNA ratio.

It is quite popular https://www.nature.com/articles/s41586-018-0414-6,  
however, the original pipeline is not well supported:
https://github.com/velocyto-team/velocyto.R/issues

There is a new one from kallisto team:
https://bustools.github.io/BUS_notebooks_R/velocity.html

# 1. Install R4.0 (development version) on O2
- module load gcc/6.2.0
- installed R-devel: https://www.r-bloggers.com/r-devel-in-parallel-to-regular-r-installation/  
because one of the packages wanted R4.0
- configure R with `./configure --enable-R-shlib` for rstudio
- remove conda from PATH to avoid using its libcurl
- module load boost/1.62.0
- module load hdf5/1.10.1
- installing velocyto.R: https://github.com/velocyto-team/velocyto.R/issues/86

# 2. Install velocyto.R on a laptop (Fedora 30 example)
bash:
```
sudo dnf update R
sudo dnf install boost boost-devel hdf5 hdf5-devel
git clone https://github.com/velocyto-team/velocyto.R
```
rstudio/R:
```
BiocManager::install("pcaMethods")
setwd("/where/you/cloned/velocyto.R")
devtools::install_local("velocyto.R")
```

# 3. Generate reference files
- `Rscriptdev `[01_get_velocity_files.R](https://github.com/naumenko-sa/crt/blob/master/velocity/01_get_velocity_files.R)
- output:
```
cDNA_introns.fa
cDNA_tx_to_capture.txt
introns_tx_to_capture.txt
tr2g.tsv
```

# 4. Index reference
This step takes ~1-2h and 100G or RAM:  
`sbatch `[02_kallisto_index.sh](https://github.com/naumenko-sa/crt/blob/master/velocity/02_kallisto_index.sh)

- inDrops3 support: https://github.com/BUStools/bustools/issues/4

# 5. References

https://www.kallistobus.tools/tutorials
https://github.com/satijalab/seurat-wrappers/blob/master/docs/velocity.md

# Velocity analysis in Python:
http://velocyto.org/
