---
title: Seurat Markers
description: This code is for finding Seurat markers
category: research
subcategory: scrnaseq
tags: [differential_analysis]
---

```bash
ssh -XY username@o2.hms.harvard.edu

srun --pty -p interactive -t 0-12:00 --x11 --mem 128G /bin/bash

module load gcc/6.2.0 R/3.4.1 hdf5/1.10.1

R
```

```r
library(Seurat)
library(tidyverse)

set.seed(1454944673L)
data_dir <- "data" 
seurat <- readRDS(file.path(data_dir, "seurat_tsne_all_res0.6.rds"))
```

Make sure the TSNEPlot looks as expected

```r
TSNEPlot(seurat)
```

Check markers for any particular cluster against all others

```r
cluster14_markers <- FindMarkers(object = seurat, ident.1 = 14, min.pct = 0.25)
```

Or look for markers of every cluster against all others

```r
seurat_markers <- FindAllMarkers(object = seurat, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
```

>**NOTE:** The `seurat_markers` object with be a dataframe with the row names as Ensembl IDs; however, since row names need to be unique, if a gene is a marker for more than one cluster, then Seurat will add a number to the end of the Ensembl ID. Therefore, do not use the row names as the gene identifiers. Use the `gene` column.

Save the markers for report generation

```r
saveRDS(seurat_markers, "data/seurat_markers_all_res0.6.rds")
```
