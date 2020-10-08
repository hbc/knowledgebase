## Need to go from counts in a Seurat object to the 10X format?

Recently found that some tools like Scrublet (doublet detection), require scRNA-seq counts to be in 10X format.
* barcodes.tsv
* genes.tsv
* matrix.mtx

How to do this easily?

```r
library(Seurat)
library(tidyverse)
library(rio)
library(DropletUtils)

# Read in Seurat object
seurat_stroma <- readRDS("./seurat_stroma_replicatePaper.rds")

# Output data
write10xCounts(x = seurat_stroma@assays$RNA@counts, path = "./cell_ranger_data_format_test")

```
