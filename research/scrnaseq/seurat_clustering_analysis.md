# Seurat singlecell RNA-Seq clustering analysis

This is a clustering analysis workflow to be run mostly on O2 using the output from the QC which is the `bcb_filtered` object. This workflow incorporates **Lorena's scripts** available within this same `Rscripts` folder. 

## Creating Seurat object at the end of the QC analysis

The first thing needed is to convert the `bcb_filtered` object in the QC to a Seurat object. We can do this by running Lorena's `bcb_to_seurat.R` script at the end of the QC analysis. The contents of the script are described below.

### Setting up the parameters

We need to load the `bcbioSingleCell` library and specify the appropriate organism and data directory to store output:

```r
library(bcbioSingleCell)

species <- "mus musculus" # change to appropriate species

data_dir <- "data"
```

### Define the cell cycle markers and save to file along with rowData

In the clustering analysis, we need to determine the likely phase of the cell cycle for each cell, to do this we need a list of markers for our organism output from the bcbioSingleCell package:

```r
cell_cycle_markers <- bcbioSingleCell::cellCycleMarkers[[camel(species)]]

s_genes <- cell_cycle_markers %>%
    filter(phase == "S") %>%
    pull("geneID")

g2m_genes <- cell_cycle_markers %>%
    filter(phase == "G2/M") %>%
    pull("geneID")

save(g2m_genes, s_genes, file = file.path(data_dir,"cycle.rda"))
```

Now save the rowData of the `bcb_filtered` data to file:

```r
saveRDS(rowData(bcb_filtered), file = file.path(data_dir,"rowData.rds"))
```

### Create Seurat object

To create the Seurat object we need only our `bcb_filtered` object, which contains the raw counts from the cells that have passed our quality control filtering parameters:

```r
seurat_raw <- CreateSeuratObject(raw.data = counts(bcb_filtered), 
                             meta.data = metrics(bcb_filtered))
                             
saveRDS(seurat_raw, file = file.path(data_dir,"seurat_raw.rds"))
```

Now that we have the Seurat object created, we can move on to the Seurat clustering analysis on O2.

## Setting up O2 environment to run clustering analysis

To run the clustering analysis on O2, be sure to have X11 forwarding working if you want to visualize any of the images. To do this, you may need to have XQuartz running on your local machine and log onto O2 with the terminal:

```bash
ssh -XY username@o2.hms.harvard.edu
```

Edit the your `.Renviron` file to have the following inside:

```bash
vim ~/.Renviron


R_LIBS_USER="/n/data1/cores/bcbio/R/library/3.4-bioc-release/library"

R_MAX_NUM_DLLS=150
```



Then start an interactive session with extra memory and x11:

```bash
srun --pty -p interactive -t 0-12:00 --x11 --mem 96G /bin/bash
```

After starting the interactive session, load the necessary R modules and start R:

```bash
module load gcc/6.2.0 R/3.4.1 hdf5/1.10.1

R
```

## Running the Seurat clustering analysis on O2

The next step is performing the actual clustering analysis with Seurat on O2. We are following **Lorena's `clustering_seurat.R` script** with descriptions and some additional plots added in. We are also using some descriptions from the **bcbioSingleCell clustering template**.

This workflow is adapted from the following sources:

- Satija Lab: [Seurat v2 Guided Clustering Tutorial](http://satijalab.org/seurat/pbmc3k_tutorial.html)
- Paul Hoffman: [Cell-Cycle Scoring and Regression](http://satijalab.org/seurat/cell_cycle_vignette.html)

To identify clusters, the following steps will be performed:

1. Normalization and transformation of the raw gene counts per cell to account for differences in sequencing depth.
2. Identification of high variance genes.
3. Regression of sources of unwanted variation (e.g. number of UMIs per cell, mitochondrial transcript abundance, cell cycle phase).
4. Identification of the primary sources of heterogeneity using principal component (PC) analysis and heatmaps.
5. Clustering cells based on significant PCs (metagenes).


### Setting up the R environment

Load the necessary libraries:

```r
library(Seurat)
library(tidyverse)
```

Create variable for where you store the needed data and load the cell cycle file stored in this directory:

```r
data_dir <- "data" 

load(file.path(data_dir, "cycle.rda")) 

set.seed(1454944673L)   

# Load Seurat object created                                                                                               

seurat_raw <- readRDS(file.path(data_dir, "seurat_raw.rds"))
```

### Subsetting to a single sample

Often identifying cell types is easiest for a single sample type. To subset the Seurat object, we can use the `SubsetData()` function. For example:

```r
pre_regressed_white <- SubsetData(pre_regressed_seurat, 
                                cells.use = rownames(pre_regressed_seurat@meta.data[which(pre_regressed_seurat@meta.data$interestingGroups == "white")])
```

### Normalizing counts, finding variable genes, and scaling the data

The raw counts are normalized using global-scaling normalization with the `NormalizeData()` function, which performs the following:

1. normalizes the gene expression measurements for each cell by the total expression 
2. multiplies this by a scale factor (10,000 by default)
3. log-transforms the result

```r
# Normalize counts for total cell expression and take log value                            

pre_regressed_seurat <- seurat_raw %>%
                        NormalizeData(normalization.method = "LogNormalize",
                                   scale.factor = 10000)  
```

Following normalization, the most variable genes are identified and will be used for downstream clustering analyses. The `FindVariableGenes()` function is called, which performs the following calculations:

1. calculates the average expression and dispersion for each gene
2. places these genes into bins
3. calculates a z-score for dispersion within each bin

This helps control for the relationship between variability and average expression. 

```r
# Find variable genes based on the mean-dispersion relationship based on z-score for dispersion. 

pre_regressed_seurat <-  pre_regressed_seurat %>%
                          FindVariableGenes(
                            mean.function = ExpMean,
                            dispersion.function = LogVMR,
                            do.plot = FALSE)
```

It's recommended to set parameters as to mark visual outliers on dispersion plot - default parameters are for ~2,000 variable genes.

Finally, the genes are scaled and centered using the `ScaleData()` function.

```r
# Scale and center data

pre_regressed_seurat <- pre_regressed_seurat %>%
                        ScaleData(model.use = "linear")
```

We can plot dispersion (a normalized measure of to cell-to-cell variation) as a function of average expression for each gene to identify a set of high-variance genes. To check that the dispersions behave as expected, decreasing with increasing mean, and to identify the most variable genes, we can visualize the dispersions with the `VariableGenePlot()` function.

```r
# Plot variable genes

VariableGenePlot(pre_regressed_seurat)
```

We can also check the number of variable genes:

```r
# Check number of variable genes to determine if correct parameters used  

length(x = pre_regressed_seurat@var.genes)
```

### Examining sources of variation in the data

Your single-cell dataset likely contains "uninteresting" sources of variation. This can include technical noise, batch effects, and/or uncontrolled biological variation (e.g. cell cycle). We can use PCA to identify these sources of variation, which can then be regressed out prior to further analysis.

### Cell cycle scoring

If we want to examine cell cycle variation in our data, we assign each cell a score, based on its expression of G2/M and S phase markers. These marker sets should be anticorrelated in their expression levels, and cells expressing neither are likely not cycling and in G1 phase. We assign scores in the `CellCycleScoring()` function, which stores S and G2/M scores in `seurat@meta.data`, along with the predicted classification of each cell in either G2M, S or G1 phase.

```r
# Perform cell cycle scoring

pre_regressed_seurat <- CellCycleScoring(
  pre_regressed_seurat,
  g2m.genes = g2m_genes,
  s.genes = s_genes)
```

Here we are checking to see if the cells are grouping by cell cycle. If we don't see clear grouping of the cells into `G1`, `G2M`, and `S` clusters on the PCA plot, then it is recommended that we don't regress out cell-cycle variation. When this is the case, remove `S.Score` and `G2M.Score` from the variables to regress (`vars_to_regress`) in the R Markdown YAML parameters.

```r
# Perform PCA and color by cell cycle phase

pre_regressed_seurat = RunPCA(
  pre_regressed_seurat,
  pc.genes = c(s_genes, g2m_genes),
  do.print = FALSE)

PCAPlot(pre_regressed_seurat, group.by= "Phase")
```

Now save the pre-regressed Seurat object:

```r
# Save pre-regression Seurat object

saveRDS(pre_regressed_seurat, file = file.path(data_dir, "seurat_pre_regress.rds"))
```

## Apply regression variables

Here we are regressing out variables of uninteresting variation, using the `vars.to.regress` argument in the `ScaleData()` function. When variables are defined in the `vars.to.regress` argument, [Seurat][] regresses them individually against each gene, then rescales and centers the resulting residuals.

We generally recommend minimizing the effects of variable read count depth (`nUMI`) and mitochondrial gene expression (`mitoRatio`) as a standard first-pass approach. If the differences in mitochondrial gene expression represent a biological phenomenon that may help to distinguish cell clusters, then we advise not passing in `mitoRatio` here.

When regressing out the effects of cell-cycle variation, include `S.Score` and `G2M.Score` in the `vars.to.regress` argument. Cell-cycle regression is generally recommended but should be avoided for samples containing cells undergoing differentiation.

```r
# Regress out the uninteresting sources of variation in the data

vars_to_regress <- c("nUMI", "S.Score", "G2M.Score")

seurat <- ScaleData(pre_regressed_seurat, vars.to.regress = vars_to_regress)
```

Now that regression has been applied, let's recheck to see if the cells are no longer clustering by cycle. We should now see the phase clusters superimpose.

```r
# Re-run the PCA plots and color by cell cycle phase

seurat <- RunPCA(
  seurat,
  pc.genes = c(s_genes, g2m_genes),
  do.print = FALSE)
  
PCAPlot(seurat, group.by= "Phase")
```

## Linear dimensionality reduction

Next, we perform principal component analysis (PCA) on the scaled data with `RunPCA()`. By default, the genes in `seurat@var.genes` are used as input, but can be defined using the `pc.genes` argument. `ProjectPCA()` scores each gene in the dataset (including genes not included in the PCA) based on their correlation with the calculated components. Though we don't use this further here, it can be used to identify markers that are strongly correlated with cellular heterogeneity, but may not have passed through variable gene selection.  The results of the projected PCA can be explored by setting `use.full = TRUE` for `PrintPCA()`.

```r
# Perform the scoring for all genes

seurat <- seurat %>%
  RunPCA(do.print = FALSE) %>%
  ProjectPCA(do.print = FALSE)
```

## Determine statistically significant principal components

To overcome the extensive technical noise in any single gene for scRNA-seq data, [Seurat][] clusters cells based on their PCA scores, with each PC essentially representing a "metagene" that combines information across a correlated gene set. Determining how many PCs to include downstream is therefore an important step. To accomplish this, we plot the standard deviation of each PC as an elbow plot with our `plotPCElbow()` function.

PC selection — identifying the true dimensionality of a dataset — is an important step for [Seurat][], but can be challenging/uncertain. We therefore suggest these three approaches to consider:

1. Supervised, exploring PCs to determine relevant sources of heterogeneity, and could be used in conjunction with GSEA for example.
2. Implement a statistical test based on a random null model. This can be time-consuming for large datasets, and may not return a clear PC cutoff.
3. **Heuristic approach**, using a metric that can be calculated instantly.

We're using a heuristic approach here, by calculating where the principal components start to elbow. The plots below show where we have defined the principal compoment cutoff used downstream for dimensionality reduction. This is calculated automatically as the larger value of:

1. The point where the principal components only contribute 5% of standard deviation (bottom left).
2. The point where the principal components cumulatively contribute 90% of the standard deviation (bottom right).

This methodology is also commonly used for PC covariate analysis on bulk RNA-seq samples.

```r
# Create elbow plot

PCElbowPlot(seurat)

# Determine the estimate for significant PCs

pct = seurat@dr$pca@sdev / sum(seurat@dr$pca@sdev) * 100
cum = cumsum(pct)
co1 = which(cum > 90 & pct < 5)[1]
co2 = sort(which((pct[1:length(pct)-1] - pct[2:length(pct)]) > 0.1),
           decreasing = T)[1] + 1 # last point where change of % of variation is more than 0.1%.
pcs = min(co1, co2) # change to any other number
```

## Cluster the cells

Seurat uses a graph-based clustering approach, inspired by SNN-Cliq [@Xu2015-je] and PhenoGraph [@Levine2015-hr]. This approach embeds cells in a graph structure, by default using a K-nearest neighbor (KNN) graph, with edges drawn between cells with similar gene expression patterns, and then attempt to partition this graph into highly interconnected ‘quasi-cliques’ or ‘communities’. As in PhenoGraph, [Seurat][] first constructs a KNN graph based on the euclidean distance in PCA space, and refines the edge weights between any two cells based on the shared overlap in their local neighborhoods (Jaccard distance). To cluster the cells, it then applies modularity optimization techniques [@Blondel2008-rf], to iteratively group cells together, with the goal of optimizing the standard modularity function.

The `FindClusters()` function implements the procedure, and contains a `resolution` argument that sets the "granularity" of the downstream clustering, with increased values leading to a greater number of clusters. We find that setting this parameter between `0.6`-`1.2` typically returns good results for single cell datasets of around 3K cells. Optimal resolution often increases for larger datasets. The clusters are saved in the `seurat@ident` slot.

Regarding the value of the `resolution` argument, use a value < 1 if you want to obtain fewer clusters.

```r
# Find cell clusters

seurat <- FindClusters(
  seurat,
  dims.use = 1:pcs,
  force.recalc = TRUE,
  print.output = TRUE,
  resolution = 0.8,
  save.SNN = TRUE)
```
A useful feature in [Seurat][] v2.0 is the ability to recall the parameters that were used in the latest function calls for commonly used functions. For `FindClusters()`, the authors provide the function `PrintFindClustersParams()` to print a nicely formatted formatted summary of the parameters that were chosen.

```r
PrintFindClustersParams(seurat)
```

# Run non-linear dimensional reduction

## t-SNE

[Seurat][] continues to use t-distributed stochastic neighbor embedding (t-SNE) as a powerful tool to visualize and explore these datasets. While we no longer advise clustering directly on t-SNE components, cells within the graph-based clusters determined above should co-localize on the t-SNE plot. This is because the t-SNE aims to place cells with similar local neighborhoods in high-dimensional space together in low-dimensional space. As input to the t-SNE, we suggest using the same PCs as input to the clustering analysis, although computing the t-SNE based on scaled gene expression is also supported using the `genes.use` argument.

```r
# Run the TSNE and plot

seurat <- RunTSNE(
  seurat,
  dims.use = 1:pcs,
  do.fast = TRUE)
  
TSNEPlot(object = seurat)
```

```r
# Save clustered cells

saveRDS(seurat, file = file.path(data_dir, "seurat_tsne.rds"))
```


> ***NOTE:***
> - *Use dev.off() if you want to save the figures generated.*
> - *Use the saved Seurat objects on a local computer to make report with figures.*
> - *rsync your data if you work on the cluster and local computer with the same data.*
