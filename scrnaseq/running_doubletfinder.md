# Running DoubletFinder


[DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder) is one of the most popular doublet finding methods with over 1200 citations since 2019 (as of Sept 2023).

## Preparing to run DoubletFinder

The key notes for running doubletFinder are:
-Each sample MUST be run separately
-Various parameters can be tweaked in the run (see doubletfinder website for details) but the most critical is the prior value of the percentage of doublets.
-DoubletFinder is not fast so best to run on O2 and save output as an RDS file.


## Step 1 - Generate subsets

Starting with your post-qc seurat object separate out each sample. Then make a list of these new objects and a vector of object names.

```r
sR01 <- subset(x = seurat_qc, subset = orig.ident %in%  c("R01"))
sW01 <- subset(x = seurat_qc, subset = orig.ident %in%  c("W01"))
s3N00 <- subset(x = seurat_qc, subset = orig.ident %in%  c("3N00"))
subsets = list(sR01,sW01,s3N00)
names = c('sR01','sW01',"s3N00")
```

## Step 2 - Run loop

This is the most computationally intensive step. Here we will loop through the list we created and run doublet finder.

```r
for (i in seq(1,length(subsets)) {

# SCT Transform and Run UMAP
obj <-  subsets[[i]]
obj <- SCTransform(obj)
obj <- RunPCA(obj)
obj <- RunUMAP(obj, dims = 1:10)

#Run doublet Finder
sweep.res.list_obj <- paramSweep_v3(obj, PCs = 1:10, sct = TRUE)
sweep.stats_obj <- summarizeSweep(sweep.res.list_obj, GT = FALSE)
bcmvn_obj <- find.pK(sweep.stats_obj)
nExp_poi <- round(0.09*nrow(obj@meta.data))  ## Assuming 9% doublet formation rate can be changed.
obj <- doubletFinder_v3(obj, PCs = 1:10, pN = 0.25, pK = 0.1, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

# Rename output columns from doublet finder to be consistent across samples for easy merging
colnames(obj@meta.data)[c(22,23)] <- c("pANN","doublet_class") ## change coordinates based on your own metadata size
assign(paste0(names[i]), obj) 
}
```

## Step 3 - Merge doubletfinder output and save

After doubletfinder is run merge to output to a single seuarat object and save that as an RDS.

```r
seurat_doublet <- merge(x = subsets[[1]], 
              y = subsets[2:length(subsets)])
saveRDS(seurat_doublet, file = "seurat_postQC_doubletFinder.rds")
```

## Step 4 (optional) add doublet info to a pre-existing seurat object for plotting

If you have gone ahead and run most of the seurat pipeline before running DoubletFinder you can add the doublet information to any object for plotting on a UMAP

```r
doublet_info <- seurat_doublet@meta.data$doublet_class
names(doublet_info) <- colnames(x = seurat_doublet)
seurat_norm <- AddMetaData(seurat_norm, metadata=doublet_info, col.name="doublet")
```

## Step 5 remove doublets

You can remove doublets from any seurat object that has the doublet info.

```r
seurat_qc_nodub <- subset(x = seurat_doublet, subset = doublet == "Singlet")
saveRDS(seurat_qc_nodub, file = "seurat_qc_nodoublets.rds")
```

