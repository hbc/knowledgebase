#Billingsley 3_7_25

#Creating a Seurat Object with Spatial data taken from Baysor


library(Seurat)
library(tidyverse)
library(Matrix)
library(RCurl)
library(knitr)
library(SeuratDisk)
library(patchwork)
library(reshape2)
library(data.table)  # Fast data loading


# Load the transcript-level data from Baysor output (segmentation.csv)

Test13<- fread("~/Desktop/Test13/segmentation.csv")

length(rownames(Test13))#3e6


#there are a number of metrics included for each transcript in the segmentation.csv file. "Confidence" is a probability that the transcript is real and not noise. "Assignment_confidence" is a probability that a transcript is unambiguously associated with a specific cell. (It's possible that a high "confidence" transcript might not have been able to have been assigned to one particular cell.) "Is_noise" is a call about whether the transcript is noise. 

#A python script written by 10x to convert Baysor output into mtx files for seurat, only used the "assignment_confidence" metric for filtering, and at a threshold of 0.9.  I think it is reasonable to use all of these and adjust as needed, thinking about sensitivity and specificity and keeping an eye on the number of transcripts being filtered. You may have to adjust these iteratively depending on downstream results.

#Here I start by filtering on is_noise

noiseF<-Test13[which(Test13$is_noise == FALSE),]

noiseT<-Test13[which(Test13$is_noise == TRUE),]

table(Test13$is_noise)

FALSE    TRUE 
2881481  251812 


hist(noiseF$confidence)# confidence assigns a probability that a transcript is real

hist(noiseT$confidence)#NoiseT has most transcripts very close to confidence 0

Test13 <-Test13 %>% filter(is_noise == FALSE)

table(Test13$is_noise)

hist(Test13$assignment_confidence)# $assignment_confidence assigns a probability that a transcript is associated with a particular called cell.

length(which(Test13$assignment_confidence >= 0.5))

length(which(Test13$assignment_confidence >= 0.9))

length(which(Test13$assignment_confidence >= 0.99))# would remove half the transcripts


Test13 <- Test13 %>% filter(assignment_confidence >= 0.9)

#length(rownames(Test13))#2.1e6


# $confidence is the likelihood that a transcript is "real"

hist(Test13$confidence)

#length(which(Test13$confidence >= 0.5))#1901953

Test13 <- Test13 %>% filter(confidence >= 0.5)


length(rownames(Test13))#1901953 


# in the $cell column, there are blanks indicating the transcript was not assigned to a cell, i.e. is noise. Here we check if there are remaining transcripts with no cell assignment, ie blank values.

which(Test13$cell == "")#no blanks

#But if so

#Test13$cell[Test13$cell == ""] <- NA

#Test13 <- Test13 %>% filter(!is.na(cell))

#Here we create a cell by transcript count matrix

gene_counts <- Test13 %>%
  count(gene, cell) %>%
  pivot_wider(names_from = cell, values_from = n, values_fill = list(n = 0))


gene_matrix <- as.matrix(gene_counts[,-1])  # Remove gene column for matrix conversion


rownames(gene_matrix) <- gene_counts$gene

str(gene_matrix)#30249 \

sum(gene_matrix)#1901953


#Remove control genes. (You can remove them at this step or later, and consider removing high expressing "nuisance" genes like MALAT1, and maybe rpls etc. I'll leave these in for now)

CtrlTranscripts<-grep("^Sys|^Neg",rownames(gene_matrix))

gene_matrix<-gene_matrix[-CtrlTranscripts,]



#Create the Seurat object


SO_t13 <- CreateSeuratObject(counts = gene_matrix)



print(sum(GetAssayData(SO_t13, slot = "counts")))# 1850745



#read in the segmentation_stats_csv file. This includes x,y coordinate data for cell centroids, and cell metrics we can use for filtering, like cell area, and cell transcript density. We can also filter on metrics like n_Count, n_Feature, Novelty.


scs<-fread(file="~/Desktop/Test13/segmentation_cell_stats.csv")

length(rownames(scs))

range(na.omit(scs$area))#0.02 : 75760 looks good, see discussion on area below.


colnames(scs)

sum(is.na(scs$area))#26


#there are 26 "cells" with NaN values for some of the metrics, like area. 

sum(is.na(scs$density))#26

which(is.na(scs$area))


scs <- scs %>%
  filter(cell %in% colnames(SO_t13))# keep only cells called in transcript annotation

#This also filtered out the NaN rows

sum(is.na(scs$density))#0 the NaN are gone


rownames(scs)<-scs$cell


SO_t13 <- AddMetaData(SO_t13, metadata = scs)


#Filtering out poor quality cells

#Filter by *cell area*, this can be an important filtering step


hist(SO_t13$area, breaks=1000)#most segments are below 5000 px^2

range(SO_t13$density)

hist(SO_t13$density, breaks=1000, ylim=c(0,200))

range(SO_t13$density)#.00248 : 3.85


length(colnames(SO_t13)[which(SO_t13$density > 2)])#6



#thinking about cell area

#My AtoMx output area range for these data  is 457 : 75800

#These are px^2 values.

#CosMx in uses a px to uM scale of 0.12028 uM to px, so 0.12028^2  uM^2 to px^2 

#so at the low range

457 * 0.12028^2 #== 6.6uM^2

#and find the radius and diameter, area = 3.14 r^2

6.6 / 3.14 #== 2.1 uM

#sqrt(2.1) == 1.44 radius == 2.88 diameter cell, pretty small



#at the high range 

75800 * 0.12028^2 #== 1096.62 uM^2

1096.62 / 3.14 #== 349 uM

sqrt(349)#  == 18.6 radius == 37 uM

#seems about right, maybe a little small . keratinocytes can be larger cells and can be 15 to 50 microns diameter


#My data here from Baysor are similar. Filter out the low area cells but leave the high end intact

range(SO_t13$area)#1.04 to 75760


#use a 3uM diameter threshold?


3.14 * 1.5^2 #== 7.065 uM^2

7.065 / 0.12028^2 #== 488 cutoff. There does appear to be a modest inflection point around there at 500

VlnPlot(SO_t13, features="area", pt.size=0) +
  ylim(0,1000)

VlnPlot(SO_t13, features="density", pt.size=0) +
  ylim(0,0.1)#probably just remove the high outliers


#What about $avg_confidence and $avg_assignment_confidence in the metadata coming from the segmentation_cell_stats.csv file? These averages were assigned to each called cell in Baysor *before* I manually filtered transcripts above, so are not accurate any longer per cell in the Seurat Object, but you could conceptually use these anyway for additional filtering.


#some additional plots

countMat<-data.frame(GetAssayData(SO_t13))

MeanVec<-apply(countMat, 1, mean)#mean counts per feature

SrtMeanVec<-sort(MeanVec, decreasing=T)

head(SrtMeanVec)#"MALAT1", "MHC I", "B2M" these are commonly seen as high expressors,  also KRT1 and KRT10 in my data. Think about removing some of the "nuisance" genes before running FAM. (Or before PCA?)

cmc<-data.frame(GetAssayData(SO_t13, assay="RNA", layer="counts"))

#features are in rows, let's count number of positives per row, i.e. positive cells pre feature. Do we need to remove any poorly expressed features?

RS<-rowSums(cmc > 0)

range(RS)# "No feature has fewer than 197 positive cells in my data" 

CS<-colSums(cmc > 0)

range(CS)#number of positive ("on-scale) features per cell. 

length(which(CS < 4))# we will remove these low feature cells

hist(MeanVec, breaks=1000, main="mean counts per feature")#there are some outlying high-count features as noted, MALAT1, B2M etc 

plot(sort(CS), log="y", main="onscale features per cell")#Im gonna remove the cells with fewer than 4 features, (some people remove cells with fewer than 10 onscale features...)


VlnPlot(SO_t13, features = "density")

VlnPlot(SO_t13, features="area", pt.size=.1) + 
  ylim(0,2000)

SO_t13$novelty <- log10(SO_t13@meta.data$nFeature_RNA) / log10(SO_t13@meta.data$nCount_RNA)

VlnPlot(object = SO_t13, features = c("nCount_RNA"), ncol = 1, pt.size = 0.05) +
  xlab("") +
  ylab("Number of counts per cell")#

VlnPlot(object = SO_t13, features = c("nCount_RNA"), ncol = 1, pt.size = 0.05) +
  ylim(0,20)
xlab("") +
  ylab("Number of counts per cell")#some zero count cells



VlnPlot(object = SO_t13, features = c("nFeature_RNA"), ncol = 1, pt.size = 0.05) +
  xlab("") +
  ylab("Number of genes detected")

VlnPlot(object = SO_t13, features = c("nFeature_RNA"), ncol = 1, pt.size = 0.05) +
  xlab("") +
  ylab("Number of genes detected") +
  ylim(0,10) #some 0 feature cells

VlnPlot(object = SO_t13, features = c("novelty"), ncol = 1, pt.size = 0.05) +
  xlab("") +
  ylab("Number of genes detected")

SO_t13@meta.data %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=novelty)) + geom_point() + 
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  xlab("nCount") + ylab("nFeature") +
  facet_wrap(~orig.ident)


length(which(CS < 4))#5219 removing low-feature cells

BadCellIndex<-which(CS < 4)


SubS1<-subset(SO_t13, cells= BadCellIndex, inv=T)

#filtering on other criteria


fDatBay13 <-subset(SubS1, subset = nFeature_RNA > 3 & novelty > 0.7 & nCount_RNA < 1500 & nCount_RNA > 5 & density <= 1 & area > 500 )#21885

print(sum(GetAssayData(fDatBay13, slot = "counts")))#1614287 




fDatBay13<-SCTransform(fDatBay13)

fDatBay13<-RunPCA(fDatBay13, npcs = 30, features = rownames(fDatBay13))

Idents(fDatBay13)<-"cluster"#these are the clusters called by Baysor

DimPlot(fDatBay13, reduction="pca", label=T)

fDatBay13<-RunUMAP(fDatBay13, dims = 1:30, return.model=T)

DimPlot(fDatBay13, reduction="umap", group.by="cluster") +
  ggtitle("Baysor clusters")#pretty squirrely, try other umap dims

#Try a number of different umap dims, find one that shows good separation of clusters without being overly-fit and ragged
#can use FeaturePlot of known markers to help with this


FeaturePlot(fDatBay13, features="KRT1", order=T)
FeaturePlot(fDatBay13, features="KRT14", order=T)
FeaturePlot(fDatBay13, features="S100A8", order=T)
FeaturePlot(fDatBay13, features="ACTA2", order=T)
FeaturePlot(fDatBay13, features="CD8A", order=T)
FeaturePlot(fDatBay13, features="CD14", order=T)


t13u10<-RunUMAP(fDatBay13, dims = 1:10, return.model=T)

DimPlot(t13u10, reduction="umap", group.by="cluster") +
  ggtitle("Baysor clusters")

FeaturePlot(t13u10, features="KRT1", order=T)
FeaturePlot(t13u10, features="KRT14", order=T)
FeaturePlot(t13u10, features="S100A8", order=T)
FeaturePlot(t13u10, features="ACTA2", order=T)
FeaturePlot(t13u10, features="CD8A", order=T)
FeaturePlot(t13u10, features="CD14", order=T)

t13u12<-RunUMAP(fDatBay13, dims = 1:12, return.model=T)

DimPlot(t13u12, group.by="cluster")

FeaturePlot(t13u12, features="KRT1", order=T)
FeaturePlot(t13u12, features="KRT14", order=T)
FeaturePlot(t13u12, features="S100A8", order=T)
FeaturePlot(t13u12, features="ACTA2", order=T)
FeaturePlot(t13u12, features="CD8A", order=T)
FeaturePlot(t13u12, features="CCL5", order=T)
FeaturePlot(t13u12, features="CD14", order=T)
FeaturePlot(t13u12, features="IGHG1", order=T)
FeaturePlot(t13u12, features="IGKC", order=T)
FeaturePlot(t13u12, features="TPSAB1/B2", order=T)
FeaturePlot(t13u12, features="CD74", order=T)
FeaturePlot(t13u12, features="LYZ", order=T)
FeaturePlot(t13u12, features="COL1A1", order=T)


#I like dims=12

fDatBay13<-RunUMAP(fDatBay13, dims = 1:12, return.model=T)

FeaturePlot(fDatBay13, reduction="umap", features="avg_confidence")

FeaturePlot(fDatBay13, reduction="umap", features="area")

FeaturePlot(fDatBay13, reduction="umap", features="density")

DimPlot(fDatBay13, reduction="umap", group.by="cluster", label=T)#Baysor clustering found the major cell types including basal and superficial ketatinocytes c7 and c9, pericytes are c6, Bcellc c1, CD8 c3, Macros C5, fibroblasts c2

fDatBay13 <- FindNeighbors(object = fDatBay13, reduction = "umap", dims = 1:2)

fDatBay13 <- FindClusters(object = fDatBay13, algorithm = 1, resolution = c(0.01, 0.025, 0.05, 0.075, 0.1, 0.125, 0.150, 0.175, .2, .25, .3), random.seed = 42)


for (res in c(0.01, 0.025, 0.05, 0.075, 0.1, 0.125, 0.150, 0.175, .2, .25, .3)){
  cat("\n")
  cat("","Cluster resolution: ",res,"\n")
  Idents(object = fDatBay13) <- paste0("SCT_snn_res.", res)
  p <- DimPlot(fDatBay13,
               reduction = "umap",
               label = TRUE) +
    ggtitle(paste0("SCT_snn.", res))
  print(p)
  cat("\n")
}

#0.125 has 15, go with this for now, but will need a higher res to get clean pericytes.


DimPlot(fDatBay13, reduction="umap", group.by="SCT_snn_res.0.125", label=T)


#FindAllMarkers

#remove nuisance genes

GrepList<-grep(rownames(fDatBay13), pattern="^RPS|^RPL|^MALAT1")

FAMt13c15<-FindAllMarkers(fDatBay13, only.pos=T, features = rownames(fDatBay13)[-GrepList])

write.table(FAMt13c15, file="~/Desktop/FAMt13c15.txt", sep="\t", col.names=NA)

#add cell coordinates for Spatial plots

centroids <- data.frame(x = fDatBay13@meta.data[["x"]], y = fDatBay13@meta.data[["y"]], cell = colnames(fDatBay13))

cents <- CreateCentroids(centroids)

coords <- CreateFOV(coords = list("centroids" = cents), type = "centroids")

fDatBay13[["FOV"]]<-coords



Idents(fDatBay13)<-"test"

ImageDimPlot(fDatBay13, fov="FOV")

Idents(fDatBay13)<-"SCT_snn_res.0.125"

ImageDimPlot(fDatBay13, fov="FOV")

Idents(fDatBay13)<-"cluster"

ImageDimPlot(fDatBay13, fov="FOV")

Idents(fDatBay13)<-"SCT_snn_res.0.125"


ImageDimPlot(fDatBay13, fov="FOV", cells=colnames(fDatBay13)[which(fDatBay13$SCT_snn_res.0.125 == 0 | fDatBay13$SCT_snn_res.0.125 == 2 )], cols=c('red', "green")) 

ImageDimPlot(fDatBay13, fov="FOV", cells=colnames(fDatBay13)[which(fDatBay13$SCT_snn_res.0.125 == 13 | fDatBay13$SCT_snn_res.0.125 == 10)], cols=c('red', "green")) + ggtitle("Baysor")


Idents(fDatBay13)<-"cluster"

DimPlot(fDatBay13, reduction="umap", label=T)

Pericytes are in Baysor cluster 6 and knn/snn cluster 5

FeaturePlot(fDatBay13, features="ACTA2", order=T)

DimPlot(fDatBay13, reduction="umap", label=T)


Idents(fDatBay13)<-"cluster"


ImageDimPlot(fDatBay13, fov="FOV", cells=colnames(fDatBay13)[which(fDatBay13$cluster == 6)], cols="red") + ggtitle("Baysor c5")

Idents(fDatBay13)<-"SCT_snn_res.0.125"

ImageDimPlot(fDatBay13, fov="FOV", cells=colnames(fDatBay13)[which(fDatBay13$SCT_snn_res.0.125 == 5)], cols="red") + ggtitle("snn/knn c5")

#they are almost identical.















