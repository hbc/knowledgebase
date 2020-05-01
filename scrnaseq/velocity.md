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

# 2. Install velocyto.R with R3.6.3 (Fedora 30 example)
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

# 5. Split reads by sample with barcode_splitter

- merge reads from multiple flowcells first
- https://pypi.org/project/barcode-splitter/
```
barcode_splitter --bcfile samples.tsv Undetermined_S0_L001_R1.fastq Undetermined_S0_L001_R2.fastq Undetermined_S0_L001_R3.fastq Undetermined_S0_L001_R4.fastq --idxread 3 --suffix .fq
```

kallisto bus counting procedure works on per sample basis, so we need to split samples to separate fastq files, and merge samples across lanes.

- [split_barcodes.sh](https://github.com/naumenko-sa/crt/blob/master/velocity/03_split_barcodes.sh)

# 6. Count spliced and unspliced transcripts
- [kallisto_count](https://github.com/naumenko-sa/crt/blob/master/velocity/04_kallisto_count.sh)
- output:
```
spliced.barcodes.txt
spliced.genes.txt
spliced.mtx
unspliced.barcodes.txt
unspliced.genes.txt
unspliced.mtx
```

# 7. Create Seurat objects for every sample
- [create_seurat_sample.Rmd](https://github.com/naumenko-sa/crt/blob/master/velocity/05.create_seurat_sample.Rmd)
- also removes empty droplets

# 8. Merge seurat objects
- [merge_seurats](https://github.com/naumenko-sa/crt/blob/master/velocity/06.merge_seurats.Rmd)

# 9. Velocity analysis
- [velocity_analysis](https://github.com/naumenko-sa/crt/blob/master/velocity/07.velocity_analysis.Rmd)

# 10. Plot velocity picture
- [plot_velocity](https://github.com/naumenko-sa/crt/blob/master/velocity/08.plot_velocity.Rmd)

# 11. Repeat marker analysis
- [velocity_markers](https://github.com/naumenko-sa/crt/blob/master/velocity/09.velocity_markers.Rmd)

# 11. References
- https://www.kallistobus.tools/tutorials
- https://github.com/satijalab/seurat-wrappers/blob/master/docs/velocity.md
- [preprocessing influences velociy analysis](https://www.biorxiv.org/content/10.1101/2020.03.13.990069v1)

# Velocity analysis in Python:
- http://velocyto.org/
- https://github.com/pachterlab/MBGBLHGP_2019/blob/master/Figure_3_Supplementary_Figure_8/supp_figure_9.ipynb
