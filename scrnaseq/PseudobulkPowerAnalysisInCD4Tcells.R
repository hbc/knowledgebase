#Billingsley 1_30_26

# This analysis is to perform a sample size power analysis for pseudobulked sc-RNA-Seq data using a pilot dataset.

# The study will be done on pbmc from human subjects with food allergy. And we want to be able to do differential expression analyses on small subsets of cells like Tregs, Tfh13 and Tfr cells. These rare CD4+ T cell subsets can make up from 1 to 4% of the total cell population.

# For pilot data we use a single cell human pbmc dataset (GSE288147) from patients with asthma and control donors. This is a fairly large dataset having 13 samples. We initially qc'd all 13 samples to look for rare cell types, but we ultimately used a subset of these to use as pilot data, 3 donors with severe disease and 3 normal controls. 

# Because it can be challenging to find these rare CD4+ cell subsets, to simulate these cell types for the power analysis for a grant submission we will collect CD4+ T cells from each donor of the number calculated to be 4% of the total pbmc count for each.

# The initial part of this analysis details the processing of the sc data. 

# The pseudobulk sample size power analysis using DESEq2 and RNASeqPower begins approx line 429


library(Seurat)
library(tidyverse)
library(Matrix)
library(RCurl)
library(knitr)
library(reticulate)
library(ensembldb)
library(AnnotationHub)
library(SeuratDisk)
library(patchwork)
library(DESeq2)
library(harmony)
library(RNASeqPower)


#Process the data in Seurat

base_dir <- "~/Desktop/GSE288147_RAW/"
LF <- list.files(base_dir)


# Prepare a list to hold Seurat objects
SOlist <- vector("list", length = length(LF))

names(SOlist) <- LF

# Loop through each directory
for (i in seq_along(LF)) {
  
  
  # Print progress
  message("Processing sample ", i, " of ", length(LF), ": ", LF[i])
  
  
  # Construct the path
  dpath <- file.path(base_dir, LF[i])
  
  # Read the 10X count data
  counts <- Read10X(data.dir = dpath)
  
  # Create a Seurat object
  so <- CreateSeuratObject(counts = counts, project = LF[i])
  
  # Add a sample ID column in metadata
  so$sample <- LF[i]
  
  # Store in the list
  SOlist[[i]] <- so
}



all_combined <- merge(
  x = SOlist[[1]],
  y = SOlist[-1],
  add.cell.ids = names(SOlist)
)


AC<-all_combined


#QC and removing poor quality droplets

AC$mitoRatio <- PercentageFeatureSet(object = AC, pattern = "^MT-")


log10GenesPerUMI <- log10(AC@meta.data$nFeature_RNA) / log10(AC@meta.data$nCount_RNA)

AC$novelty <- log10GenesPerUMI

AC@meta.data %>%   
  ggplot() + 
  geom_bar(aes(x=orig.ident, fill=orig.ident)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust=0.5)) +
  ggtitle("Number of Cells") +
  xlab("")


table(AC$orig.ident)

Idents(AC)<-"orig.ident"

VlnPlot(object = AC, features = c("nCount_RNA"), ncol = 1, pt.size = 0.0) +
  xlab("") +
  ylab("Number of umi per cell")



VlnPlot(object = AC, features = c("nCount_RNA"), ncol = 1, pt.size = 0.0) +
  xlab("") +
  ylab("Number of umi per cell") + ylim(0,20000) #use maybe 15000

VlnPlot(object = AC, features = c("nCount_RNA"), ncol = 1, pt.size = 0.0) +
  xlab("") +
  ylab("Number of umi per cell") + ylim(0,1000)

VlnPlot(object = AC, features = c("nFeature_RNA"), ncol = 1, pt.size = 0.0) +
  xlab("") +
  ylab("Number of umi per cell") + ylim(0,20000)

VlnPlot(object = AC, features = c("nFeature_RNA"), ncol = 1, pt.size = 0.0) +
  xlab("") +
  ylab("Number of umi per cell") + ylim(0,1000)#some droplets with pretty low feature numbers, like peaks around100, use 500

VlnPlot(object = AC, features = c("mitoRatio"), ncol = 1, pt.size = 0.0) +
  xlab("") +
  ylab("Number of umi per cell") + ylim(0,30)#


VlnPlot(object = AC, features = c("novelty"), ncol = 1, pt.size = 0.0) +
  xlab("") 


AC@meta.data %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=novelty)) + 
  geom_point() + 
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  xlab("nCount_RNA") + ylab("nFeature_RNA") +
  facet_wrap(~orig.ident)


#There are *many* poor quality cells


#Filter poor quality cells.


fD<-subset(AC, subset = nCount_RNA < 15000  & nFeature_RNA > 500 & novelty > 0.8 & mitoRatio < 15 ) 

fD@meta.data %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=novelty)) + 
  geom_point() + 
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  xlab("nCount_RNA") + ylab("nFeature_RNA") +
  facet_wrap(~orig.ident)

length(colnames(AC))#133k


length(colnames(fD))#105k

AC@meta.data %>%   
  ggplot() + 
  geom_bar(aes(x=orig.ident, fill=orig.ident)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust=0.5)) +
  ggtitle("Number of Cells") +
  xlab("") +
  ylim(0,22000)

table(AC$orig.ident)

fD@meta.data %>%   
  ggplot() + 
  geom_bar(aes(x=orig.ident, fill=orig.ident)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust=0.5)) +
  ggtitle("Number of Cells") +
  xlab("") +
  ylim(0,22000)

table(fD$orig.ident)


saveRDS(fD, file="~/Desktop/fD_1_26_26.rds")



#Normalize the data


fD<-SCTransform(fD)

fD<-RunPCA(fD)

fD<-RunUMAP(fD, dims=c(1:20), return.model=T)

DimPlot(fD, reduction="umap")


#Data need to be Harmony integrated for cell annotation


fD<-RunHarmony(fD, "orig.ident", plot_convergence = TRUE, assay.use="SCT", reduction="pca", reduction.save="harmonyRNA")

fD<-RunUMAP(fD, reduction = "harmonyRNA", dims=1:10, reduction.name ="harmonyUMAP", reduction.key = "harmonyUMAP_")

#Identify some basic pbmc cell types


DimPlot(fD, reduction="harmonyUMAP")
FeaturePlot(fD, reduction="harmonyUMAP", features="CD3G")

DimPlot(fD, reduction="harmonyUMAP")
FeaturePlot(fD, reduction="harmonyUMAP", features="CD14")

DimPlot(fD, reduction="harmonyUMAP")
FeaturePlot(fD, reduction="harmonyUMAP", features="FCGR3A")


#perform basic snn/knn cluster identification

fD <- FindNeighbors(object = fD, reduction = "pca", dims=1:20)

fD<-FindClusters(object = fD, algorithm = 1, resolution = c(0.01, .025, 0.035, .05, .1, .15, .2), random.seed = 42)

for (res in c(0.01, 0.025, 0.035, 0.05, .1, .15, .2)){
  cat("\n")
  cat("","Cluster resolution: ",res,"\n")
  Idents(object = fD) <- paste0("SCT_snn_res.", res)
  p <- DimPlot(fD,
               reduction = "harmonyUMAP",
               label = TRUE) +
    ggtitle(paste0("SCT_snn.", res))
  print(p)
  cat("\n")
}

#use 0.01

Idents(fD)<-"SCT_snn_res.0.01"

fD<-FindClusters(object = fD, algorithm = 1, resolution = 0.01, random.seed = 42)


#perform basic differential expression analyses to help id pbmc cell populations

DefaultAssay(fD)<-"SCT"


GrepList<-grep(rownames(fD), pattern="^RPS|^RPL|MALAT1")#remove nuisance genes

length(GrepList)#101 


#collect 10% of each cluster to use as a reduced dataset that will run faster, this is fine for annotating the basic cell clusters, and speed things up

cells_to_keep <- fD@meta.data %>%
  rownames_to_column("cell") %>%
  group_by(SCT_snn_res.0.01) %>%
  # sample 10% of cells from each cluster
  slice_sample(prop = 0.10) %>%
  pull(cell)

seurat_subset <- subset(fD, cells = cells_to_keep)



FAM<-FindAllMarkers(seurat_subset[-GrepList], only.pos=T, recorrect_umi = FALSE)

top100_by_padj <- FAM %>%
  group_by(cluster) %>%
  slice_min(order_by = p_val_adj, n = 100, with_ties = FALSE) %>%
  ungroup()

write.table(top100_by_padj, "~/Desktop/FAM_1_26.txt", sep="\t", col.names=NA)


FeaturePlot(fD, reduction="harmonyUMAP", features="FOXP3")

#FOXP3 + cells dont aggregate well in the harmonized object

#we decided on using small numbers of CD4+ T cells as a surrogate for Tregs.

FeaturePlot(fD, reduction="harmonyUMAP", features="CD8A")


#Subset the T cells, then the CD4+ T cells.


p1<-DimPlot(fD, reduction="harmonyUMAP")


#Grab the T cells

Tcells<-CellSelector(p1)

Tsub<-subset(fD, cells=Tcells)


#Recluster the T cells

Tsub<-RunPCA(Tsub)



Tsub<-RunHarmony(Tsub, "orig.ident", plot_convergence = TRUE, assay.use="SCT", reduction="pca", reduction.save="harmonyRNA")

Tsub<-RunUMAP(Tsub, reduction = "harmonyRNA", dims=1:20, reduction.name ="harmonyUMAP", reduction.key = "harmonyUMAP_")

DimPlot(Tsub, reduction="harmonyUMAP")

FeaturePlot(Tsub, reduction="harmonyUMAP", features="FOXP3")

FeaturePlot(Tsub, reduction="harmonyUMAP", features="CD8A")

FeaturePlot(Tsub, reduction="harmonyUMAP", features="CD4")

FeaturePlot(Tsub, reduction="harmonyUMAP", features="NCAM1")

FeaturePlot(Tsub, reduction="harmonyUMAP", features="FCGR3A")

saveRDS(Tsub, file="~/Desktop/Tsub_1_26.rds")


table(Tsub$orig.ident)

Idents(Tsub)<-"orig.ident"

#Here we decide to collect a subset of samples for the power analysis. We initially tried using POWSC which was having memory issues with larger cell numbers...ultimately we decided to go with pseudobulk approach with RNASeqPower, as POWSC apparently doesnt do sample size estimates, just cluster cell number estimates, but I liked the idea of a reduced dataset for that as well.

#collect three "severe" and three "control" samples

Tsub6<-subset(Tsub, idents=c("B02", "B03", "B04", "C02", "C03", "C04"))


#Here, we will use snn/knn clusters to identify the CD4+ T cells from the remaining T cells, because a square lasso using CellSelector wont work with the shape of the umap

Tsub6<-FindNeighbors(object=Tsub6, reduction="harmonyUMAP", dims=1:2)## Typically when finding neighbors, you use pca reduction to use maximum dimeensions, but that can leave fuzzy clusters. Here we are just using clustering to identify clean "gross" clusters to get the CD4+ cells so we use the harmonyUMAP reduction.


Tsub6<-FindClusters(object = Tsub6, algorithm = 1, resolution = c(0.01, .025, 0.035), random.seed = 42)

for (res in c(0.01, 0.025, 0.035)){
  cat("\n")
  cat("","Cluster resolution: ",res,"\n")
  Idents(object = Tsub6) <- paste0("SCT_snn_res.", res)
  p <- DimPlot(Tsub6,
               reduction = "harmonyUMAP",
               label = TRUE) +
    ggtitle(paste0("SCT_snn.", res))
  print(p)
  cat("\n")
}

#Grab cluster 0

Idents(Tsub6)<-"SCT_snn_res.0.01"

c0<-subset(Tsub6, idents="0")


saveRDS(c0, file="~/Desktop/c0_1_26.rds")


#c0<-readRDS(file="~/Desktop/c0_1_26.rds")




#We initially looked for Tregs, the clients celltype of interest, but these were not easily identified in the pilot data. In the interest of time we decided to use CD4+ T cells as a surrogate cell population (Tregs are CD4+ T cells). But we wanted to use a similar single cell number input as this can influence estimates of aggregated sequencing depth. As Tregs are expected to be 1-4% of pbmc we decided to calculate 4% of each of the 6 samples cell counts and use that to collect that number of CD4+ T cells from each sample.

#Calculate 4% of total pbmc numbers (fD) for each of the 6 samples


#fD<-readRDS("~/Desktop/fD_1_26_26.rds")

sub6<-subset(fD, idents=c("B02", "B03", "B04", "C02", "C03", "C04"))

table(sub6$orig.ident)

cell_counts <- c(B02 = 3081, B03 = 11579, B04 = 10265,
                 C02 = 6743, C03 = 7562, C04 = 16056)

four_percent <- round(cell_counts * 0.04)

#Collect this number of single cells from each of these samples into a new matrix from the CD4+ sprted cells (c0)

selected_cells <- character(0)

# loop through each donor
for (donor in names(four_percent)) {
  # get all cell barcodes for that donor
  donor_cells <- WhichCells(c0, expression = orig.ident == donor)
  
  # how many to sample for this donor
  n <- four_percent[donor]
  
  # sample without replacement
  set.seed(42)  # for reproducibility
  sampled_barcodes <- sample(donor_cells, size = n, replace = FALSE)
  
  # append to the master list
  selected_cells <- c(selected_cells, sampled_barcodes)
}

length(selected_cells)

sub_selected <- subset(c0, cells = selected_cells)


length(colnames(sub_selected))#2211


# sanity Check
table(sub_selected$orig.ident)

DimPlot(sub_selected, reduction="harmonyUMAP")

DimPlot(sub_selected, reduction="umap", group.by="orig.ident")







#Pseudobulk-based approach power analysis using RNASeqPower and DESeq2.

#For an edgR approach see here:

#https://www.dropbox.com/scl/fo/qj9o9t2fnnfyofqw4h6x9/AJXag0767hoOiN1_x-8Gg3I?rlkey=zcyoigteqvrpxi2hj3hpae5rs&dl=0

#For this we need to calculate input estimates of data read depth, and dispersion.


#aggregate by donor the singlecell counts to pseudobulk counts

pseudobulk <- AggregateExpression(
  object     = sub_selected,
  assays     = "RNA",        # which assay to aggregate on
  group.by   = c("orig.ident"),      # sample/replicate grouping
  return.seurat = TRUE       # return a Seurat object
)

str(pseudobulk)#6 x 36k seuratobject


#get the metadata back as AggregateExpression loses that.

#add new disease factor to object, "severe" or "control"

groupFac<-factor(sub_selected$orig.ident)

levels(groupFac)<-c("sev", "sev", "sev", "ctrl", "ctrl", "ctrl" )

table(groupFac)

sub_selected$group<-groupFac


orig_to_group <- unique(sub_selected@meta.data[, c("orig.ident", "group")])


pb_meta <- unique(orig_to_group[, c("orig.ident", "group")])

colnames(pb_meta) <- c("sample", "condition")


#get the counts matrix

pb_counts <- GetAssayData(
  object = pseudobulk,
  assay  = "RNA",
  layer  = "counts"
)

str(pb_counts)#matrix 36k by 6

# here we filter out poor quality transcripts. This might not be necessary for running DE in DESeq2, but for generating dispersion estimates for power analysis it can make a difference by removing noise.


keep <- rowSums(pb_counts >= 10) >= 3 #this keeps genes with at least counts of 10 in 3 of the 6 samples

table(keep)#9901 transcripts out of 36k

pb_counts_filtered <- pb_counts[keep, ] # 6 by 9k

rownames(pb_meta)<-colnames(pb_counts_filtered)

#run DESeq2 to get the dispersions


dds_filtered_4pc <- DESeqDataSetFromMatrix(
  countData = pb_counts_filtered,
  colData   = pb_meta,
  design    = ~ condition
)


dds_filtered_4pc <- DESeq(dds_filtered_4pc)

res_filtered_4pc <- results(dds_filtered_4pc)

summary(res_filtered_4pc)


#summarize variability...

#there are some decisions to make here, how to calculate a dispersion estimate to input to RNASeqPower.


"As the mean grows to infinity, the
square root of dispersion gives the coefficient of variation for the counts" #from the Deseq2 manual. 


#So if count means are sufficiently large, the sqrt of dispersion can be an estimate of CV.


#1. Take the median or mean of dispersions then sqrt

#as such for mean:

avg_disp_filt_4pc <- mean(dispersions(dds_filtered_4pc), na.rm = TRUE)

cv_filt_4pc_mean <- sqrt(avg_disp_filt_4pc)#0.301

#for median

med_disp_filt_4pc <- median(dispersions(dds_filtered_4pc), na.rm = TRUE)

cv_filt_4pc_med <- sqrt(med_disp_filt_4pc)#0.259



hist(as.numeric(dispersions(dds_filtered_4pc)), breaks=1000)

#The histogram is right skewed, so the mean of the dispersions will be higher than the median. This will change the final sample number calculations, using the mean will result in somewhat higher numbers of suggested samples for a given power, so mean would be a more conservative approach.





#2. an alternative more accurate approach, find the *gene-specific CVs* then take the median of those

#“The dispersion is essentially the squared coefficient of variation: Var = mu + dispersion * mu^2

#Var / mu^2 = 1/mu + dispersion
#CV² = 1/mu + dispersion” 

"When the counts are large, e.g. 100, then 1/100 is small compared to a typical dispersion value, so you have CV ~= sqrt(dispersion)" # This suggests ok to use approach 1 above if your counts are relativelt high.


#see Michael Love https://support.bioconductor.org/p/88880/

#For greater accuracy we can use CV² = 1/mu + dispersion
  
#mu <- mcols(dds)$baseMean #This gives you the mean expression of each gene

#alpha <- dispersions(dds) #Gives you the dispersion of each gene


#cv_gene <- sqrt(1/mu + alpha) # cv for each gene
  

mu<-mcols(dds_filtered_4pc)$baseMean

alpha<-dispersions(dds_filtered_4pc)


cv_gene <- sqrt(1/mu + alpha)

#then take the median value of cv_gene for all the genes

hist(as.numeric(cv_gene), breaks=1000)#it's right skewed, so use median


cv_med <- median(cv_gene[is.finite(cv_gene)], na.rm = TRUE) # 0.305 is very close to using the mean of the dispersions from approach 1.


#next we need a measure of read depth, we'll take the geometric mean of the rowmeans (mean of each gene across samples).

geomean<-exp(mean(log(rowMeans(pb_counts_filtered) + 1))) - 1

geomean#47.6



#run rnapower.



target_fc <- log2(1.5)

depth_filt    <- geomean      # typical per-gene mean

effects   <- c(1.5,2,4,8)      # targets log2 fold change

alpha    <- 0.05           # Type I error rate

target_powers <- c(0.8, 0.9)


results <- expand.grid(
  effect = effects,
  power = target_powers
)
results$samples_per_group <- NA_real_


for (i in seq_len(nrow(results))) {
  eff <- results$effect[i]
  pw  <- results$power[i]
  
  # rnapower() can solve for sample size when 'power' is given
  ss <- rnapower(
    depth  = depth_filt,
    cv     = cv_med,
    effect = eff,
    alpha  = alpha,
    power  = pw
  )
  
  # rnapower() returns a matrix/array when solving sample size
  # so take the result directly
  results$samples_per_group[i] <- as.numeric(ss)
}


results

For an edgR approach see here:

https://www.dropbox.com/scl/fo/qj9o9t2fnnfyofqw4h6x9/AJXag0767hoOiN1_x-8Gg3I?rlkey=zcyoigteqvrpxi2hj3hpae5rs&dl=0










