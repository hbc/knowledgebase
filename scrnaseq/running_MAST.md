# Running MAST

[MAST](https://github.com/RGLab/MAST) analyzes differential expression using the cell as the unit of replication rather than the sample (as is done for pseduobulk)
                                                                                                                                         
**NOTES**
                                                                                                                                         
-MAST uses a hurdle model designed for zero heavy data.
-MAST "expects" log transformed count data. 
-A [recent paper](https://doi.org/10.1038/s41467-021-21038-1) advises the use of sample id as a random factor to prevent pseudoreplication               
-Most MAST models include the total number of genes expressed in the cell.

## Where to run MAST?

MAST can be run directly in Seurat (via [FindMarkers](https://satijalab.org/seurat/reference/findmarkers) or by itself.

While MAST is easier to run with Seurat there are two big downsides:

1. Seurat does not log transform the data
2. You cannot edit the model with seurat, meaning that you cannot add sample ID or number of genes expressed.

## Running MAST (1) - from seurat to a SCA object


```r
# Seurat to SCE
sce <- as.SingleCellExperiment(seurat_obj)

# Add log counts
assay(sce, "log") = log2(counts(sce) + 1)

# Create new sce object (only 'log' count data)
sce.1 = SingleCellExperiment(assays = list(log = assay(sce, "log")))
colData(sce.1) = colData(sce)

# Change to SCA
sca = SceToSingleCellAssay(sce.1)

```

## Running MAST (2) - Filter SCA object

Here we are only filtering for genes expressed in 10% of cells but this can be altered and other filters can be added.

```r
expressed_genes <- freq(sca) > 0.1
sca_filtered <- sca[expressed_genes, ]

```

## Format SCA metadata

We add the total number of genes expressed per cell as well as setting factors as factors and scaling all continuous variables as suggested by MAST.

```r
cdr2 <- colSums(SummarizedExperiment::assay(sca_filtered)>0)
 
SummarizedExperiment::colData(sca_filtered)$ngeneson <- scale(cdr2)
SummarizedExperiment::colData(sca_filtered)$orig.ident <- factor(SummarizedExperiment::colData(sca_filtered)$orig.ident)
SummarizedExperiment::colData(sca_filtered)$Gestational_age_scaled <- scale(SummarizedExperiment::colData(sca_filtered)$Gestational_age)
```

## Run MAST

This is the most computationally instensive step and takes the longest. 
Here our model includes the number of genes epxressed (ngeneson), sample id as a random variable ((1 | orig.ident)), Gender, and Gestational age scaled.

We extract the results from our model for our factor of interest (Gestational_age_scaled)


```r
zlmCond <- suppressMessages(MAST::zlm(~ ngeneson + Gestational_age_scaled + Gender + (1 | orig.ident),  sca_filtered, method='glmer',ebayes = F,strictConvergence = FALSE))
summaryCond <- suppressMessages(MAST::summary(zlmCond,doLRT='Gestational_age_scaled'))
```

## Format Results

MAST results look quite different than DESeq2 results so we need to apply a bit of formatting to make them readable.

After formatting outputs can be written directly to csv files.
```r
summaryDt <- summaryCond$datatable

# Create reable results table for all genes tested
fcHurdle <- merge(summaryDt[contrast == "Gestational_age_scaled"
 & component == 'H', .(primerid, `Pr(>Chisq)`)], # This extracts hurdle p-values 
                  summaryDt[contrast == "Gestational_age_scaled" & component == 'logFC', 
                            .(primerid, coef, ci.hi, ci.lo)], 
                  by = 'primerid') # This extract LogFC data

fcHurdle <- stats::na.omit(as.data.frame(fcHurdle))

fcHurdle$fdr <- p.adjust(fcHurdle$`Pr(>Chisq)`, 'fdr')


# Create reable results table for significant genes
fcHurdleSig <- merge(fcHurdle[fcHurdle$fdr < .05,],
                     as.data.table(mcols(sca_filtered)), by = 'primerid')
setorder(fcHurdleSig, fdr)

```
