## Musings on T Cell Receptor (TCR) and B Cell Receptor (BCR) V,D,J and C transcripts, and their influence on clustering in single-cell RNA-Seq experiments
<br>
<br>

In vertebrates, the lymphocyte antigen receptors, T Cell Receptor (TCR) and B Cell Receptor (BCR) are encoded by V ("variable"), D ("diversity"), J ("joining"), and C ("constant") genes (or gene "segments") located in TCR and BCR loci. Numerous different versions of each segment type are found in germline DNA. (For example, in humans there are > 80 different TCR*alpha* V segments.) During lymphocyte development, through a process called "somatic recombination", or "v(d)j recombination" (in contrast to "genetic recombination" during meiosis), individual v,d,j and c segments are selected to encode a unique TCR or BCR sequence for each cell. The unselected segments are excised from the DNA. The nearly random selection and combination of individual v,d,j and c segments, promotes a vast (est. > 10<sup>15</sup> ) combinatorial diversity of potential unique antigen binding domains encoded on individual lymphocyte antigen receptors. Antigen receptor diversity generated through somatic recombination, and immune memory promoted through lymphocyte clonal expansion are the hallmarks of vertebrate adaptive immunity.


[https://en.wikipedia.org/wiki/V(D)J_recombination](https://en.wikipedia.org/wiki/V\(D\)J_recombination)

[v(d)j recombination](https://www.ncbi.nlm.nih.gov/books/NBK27140/)

Reference nomenclature for v,d,j gene loci in mouse/human: [IMGT]([http://www.imgt.org/genedb/](https://www.imgt.org/IMGTrepertoire/))

Note on the _nearly_ randomness of the v.d.j recombination process: In practice, germline v,d,j genes do not have the exact same probability of being used for a rearrangement, notably because they differ in terms of their copy number, DNA conformation/chromatin accessibility, and because the process of segment excision is iterative (if a non-functional (= non-productive) rearrangement is generated after a first round of excision favoring genes in the middle of the locus, a second round of rearrangement may occur using more distal genes before the clone is discarded). See [this review](https://pubmed.ncbi.nlm.nih.gov/29944757/) for a quantitative discussion of the matter and its consequence when assessing repertoire sharing.
<br>
<br>

**Antigen receptor transcripts and clustering**


The presence of RNA transcript counts arising from antigen receptor sequences can promote clustering in sc-RNA-Seq analysis. This can be biologically interesting, or not, depending on the experiment.

For example, there are two subtypes of TCRs, TCR*a*/*b* (alpha/beta), or TCR*g*/*d* (gamma/delta). These two different T cell subsets have distinct biological activities, so segregation of these two subsets is generally useful, and the differential expression of TRA, TRB, TRG and TRD transcripts are useful to promote subset segregation and annotation. Broadly speaking, CD4/CD8 T lymphocytes express alpha/beta TCRs while gamma/delta T lymphocytes express gamma/delta TCRs in absence of CD4/CD8 co-stimulatory protein. However, there are some exceptions: for example, a sub-population of [CD4+ Vd1+](https://pubmed.ncbi.nlm.nih.gov/25709606/) (i.e. using TCR V gene, segment delta 1) T cells. 

There are also two general subtypes of BCRs, defined by the usage of either lambda ( *l* ) or kappa ( *k* ) light chains. The biological activities, however, of these two B cell subtypes are less distinct, so clustering induced by differential *l* or *k* gene segment usage may not be interesting to a specific experiment, and in fact, can often promote significant subclustering which can interfere with identification and annotation of more biologically or experimentally interesting clustering and annotation. (Although, there may be some experimental situations where quantification of kappa and lambda B cells is of interest.)
<br>
<br>

***

### Removing antigen receptor transcripts from Seurat objects
<br>
<br>

**Option 1: By filtering the Seurat object:**
<br>
<br>

If the analyst and the client are interested in observing the effects on clustering promoted by antigen receptor genes, these are readily removed by subsetting the object.

This code displays the names of the TCR*a* chain V and J transcripts found in the Seurat object. (TCR*a* and *d* chains don't utilize D segments)

````
DefaultAssay(SeuratObject)<-"RNA"

rownames(SeuratObject)[(grep(pattern = '^TRAV|^TRAJ' , rownames(SeuratObject)))]
````
		
which will return something like this:
````
[1] "TRAV1-1"     "TRAV1-2"     "TRAV2"       "TRAV3"       "TRAV4"      
  [6] "TRAV5"       "TRAV6"       "TRAV7"       "TRAV8-1"     "TRAV9-1"    
 [11] "TRAV10"      "TRAV11"      "TRAV12-1"    "TRAV8-2"     "TRAV8-3"    
 [16] "TRAV13-1"    "TRAV12-2"    "TRAV8-4"     "TRAV8-5"     "TRAV13-2"   
 [21] "TRAV14DV4"   "TRAV9-2"     "TRAV15"      "TRAV12-3"    "TRAV8-6"    
 [26] "TRAV16"      "TRAV17"      "TRAV18"      "TRAV19"      "TRAV20"     
 [31] "TRAV21"      "TRAV22"      "TRAV23DV6"   "TRAV24"      "TRAV25"     
 [36] "TRAV26-1"    "TRAV8-7"     "TRAV27"      "TRAV28"      "TRAV29DV5"  
 [41] "TRAV30"      "TRAV31"      "TRAV32"      "TRAV33"      "TRAV26-2"   
 [46] "TRAV34"      "TRAV35"      "TRAV36DV7"   "TRAV37"      "TRAV38-1"   
 [51] "TRAV38-2DV8" "TRAV39"      "TRAV40"      "TRAV41"      "TRAJ61"     
 [56] "TRAJ60"      "TRAJ59"      "TRAJ58"      "TRAJ57"      "TRAJ56"     
 [61] "TRAJ55"      "TRAJ54"      "TRAJ53"      "TRAJ52"      "TRAJ51"     
 [66] "TRAJ50"      "TRAJ49"      "TRAJ48"      "TRAJ47"      "TRAJ46"     
 [71] "TRAJ45"      "TRAJ44"      "TRAJ43"      "TRAJ42"      "TRAJ41"     
 [76] "TRAJ40"      "TRAJ39"      "TRAJ38"      "TRAJ37"      "TRAJ36"     
 [81] "TRAJ35"      "TRAJ34"      "TRAJ33"      "TRAJ32"      "TRAJ31"     
 [86] "TRAJ30"      "TRAJ29"      "TRAJ28"      "TRAJ27"      "TRAJ26"     
 [91] "TRAJ25"      "TRAJ24"      "TRAJ23"      "TRAJ22"      "TRAJ21"     
 [96] "TRAJ20"      "TRAJ19"      "TRAJ18"      "TRAJ17"      "TRAJ16"     
[101] "TRAJ14"      "TRAJ13"      "TRAJ12"      "TRAJ11"      "TRAJ10"     
[106] "TRAJ9"       "TRAJ8"       "TRAJ7"       "TRAJ6"       "TRAJ5"      
[111] "TRAJ4"       "TRAJ3"       "TRAJ2"       "TRAJ1"  
````
<br>
<br>



As mentioned above, defining TCR *a*/*b* and TCR *g*/*d* subsets is usually interesting. These can be identified by C region transcripts.

````
rownames(SeuratObject)[(grep(pattern = '^TRAC|^TRBC|^TRGC|^TRDC' , rownames(SeuratObject)))]


[1] "TRGC2" "TRGC1" "TRBC1" "TRBC2" "TRDC"  "TRAC" 
`````


<br>
<br>

To remove all TCR*a* and TCR*d* V and J transcripts, and all TCR*b* and TCR*g* V, D, and J segments, you can subset the object. 

````
GrepX<-grep(c("^TRAV|^TRAJ|^TRBV|^TRBD|^TRBJ|^TRGV|^TRGD|^TRDJ|^TRDV|^TRDD"), rownames(SeuratObject))

NewSubsetObject<-subset(SeuratObject, features=(-GrepX))

````
<br>
<br>


Similarly, to remove all BCR transcripts use this:

````
GrepZ<-grep(c("^IGLV|^IGLJ|^IGKV|^IGKJ|^IGHV|^IGHD|^IGHJ|^IGHA|^IGHE|^IGHG|^IGHD|^IGHM|^IGLC|^IGKC"), rownames(SeuratObject))

````

The above code will accidentally collect a non-BCR gene IGHMBP2, so we'll remove that index from the index vector before subsetting.

<br>
<br>

More code

````

#There are no D regions in LC, l or k chains

#BCR Light Chain V and J transcripts


rownames(SeuratObject)[(grep(pattern = c("^IGLV|^IGLJ") , rownames(SeuratObject)))]

rownames(SeuratObject)[(grep(pattern = c("^IGKV|^IGKJ") , rownames(SeuratObject)))]



#BCR Heavy Chain, V, D, J transcripts

rownames(SeuratObject)[(grep(pattern = c("^IGHV|^IGHD|^IGHJ") , rownames(SeuratObject)))]



#BCR C regions

#Heavy Chain C regions Defining Isotypes

rownames(SeuratObject)[(grep(pattern = c("^IGHA|^IGHE|^IGHG|^IGHD|^IGHM") , rownames(SeuratObject)))]

#this includes IGHMBP2 by accident



#BCR Light Chain C regions 

rownames(SeuratObject)[(grep(pattern = c("^IGLC|IGKC") , rownames(SeuratObject)))]



#Get everything

GrepZ<-grep(c("^IGLV|^IGLJ|^IGKV|^IGKJ|^IGHV|^IGHD|^IGHJ|^IGHA|^IGHE|^IGHG|^IGHD|^IGHM|^IGLC|^IGKC"), rownames(SeuratObject))


length(GrepZ)#409

#remove the IGHMBP2 index from the GrepZ index vector

which(rownames(SeuratObject) == "IGHMBP2")#20143

GrepZ<-GrepZ[! GrepZ == 20143]



#subset

NewSubsetObject<-subset(SeuratObject, features= rownames(SeuratObject)[-GrepZ])
````
<br>
<br>

Here is a link to an .html document  demonstrating removal of BCR transcripts

[link](https://www.dropbox.com/scl/fi/lmhvjmfjw70hfa1wvf8tj/VisterraSomeCorrectedCode7_21_23.html?rlkey=ctqj9ixkmkymqnwviu5rsn4fc&dl=0)

<br>
<br>


**OPTION 2: By excluding TCR transcripts from variable features before dimensionality reduction:**
<br>
<br>

An alternative approach consists in keeping all genes in the Seurat object (thus avoiding bias in normalization due to the removal of those genes, which may be highly expressed in activated T lymphocytes and thus represent a large fraction of a cell's captured transcriptome). 

_NOTE:_ In the following, you will notice that we only removed TRAV/TRAJ and TRBV/TRBJ genes. This is because in practice TRDV/TRDJ segments are rarely sequenced (they are so short that it is difficult to map them to the genome, especially since the germline sequence is further modified by adding some non-templated nucleotides before joining the v,d and j segments).

For a single sample (no merging or integration of the data), simply add this to your standard normalization and dimensionality reduction workflow, here with SCT:

```
# Perform SCT-normalization
seurat_clean <- SCTransform(seurat_clean, verbose = FALSE,
                            variable.features.n = 3200)

# Exclude TRAV/TRAJ and TRBV/TRBJ genes from top features for PCA reduction
var_features <- VariableFeatures(seurat_clean)
delIdx <- c(grep("Trbv", var_features), grep("Trbj", var_features), 
            grep("Trav", var_features), grep("Traj", var_features))
var_features[delIdx]  # Print outcomes of `grep()` to ensure no exclusion of other genes!!
if (length(delIdx) > 0) { var_features <- var_features[-delIdx] }

# Manually force VariableFeatures() to the top 3000 features in the remaining vector
VariableFeatures(seurat_clean) <- var_features[1:3000]

# Calculate PCs
DefaultAssay(seurat_clean) <- "SCT"
seurat_clean <- RunPCA(seurat_clean, assay = "SCT", npcs = 50)

# Perform dimensionality reduction on SCT assay
DefaultAssay(seurat_clean) <- "SCT"
seurat_clean <- RunUMAP(seurat_clean, 
                        reduction = "pca", 
                        reduction.name = 'umap.rna', 
                        reduction.key = "UMAPxRNA_",
                        dims = 1:rpcs)
```


When working with multiple samples and running SCT on each independently before merge (or if applicable, integration), use the following:

```
# Perform SCT-normalization on each individual sample
seurat_ls <- SplitObject(seurat_clean, split.by = "orig.ident")

for (i in 1:length(seurat_ls)) {
  print(paste0("Processing sample ", i))
  seurat_ls[[i]] <- SCTransform(seurat_ls[[i]], verbose = TRUE,
                                variable.features.n = 3000,
                                vars.to.regress = rna_regress)
}

# Find most variable features across samples to merge
var_features <- SelectIntegrationFeatures(object.list = seurat_ls, nfeatures = 3200) 

# Exclude TRAV/TRAJ and TRBV/TRBJ genes from top features for PCA reduction
delIdx <- c(grep("TRBV", var_features), grep("TRBJ", var_features), 
            grep("TRAV", var_features), grep("TRAJ", var_features))
var_features[delIdx]
var_features <- var_features[-delIdx]
var_features <- var_features[1:3000]


# Merge SCT-normalized samples
seurat_sct <- merge(x = seurat_ls[[1]],
                    y = seurat_ls[2:length(seurat_ls)],
                    merge.data = TRUE)
DefaultAssay(seurat_sct) <- "SCT"

# Manually set variable features of merged Seurat object
VariableFeatures(seurat_sct) <- var_features

# Calculate PCs using manually set variable features
seurat_sct <- RunPCA(seurat_sct, assay = "SCT", npcs = 50)
seurat_sct <- RunUMAP(seurat_sct, 1:40)
```

For integration with Seurat CCA rather than simple sample merge, use this:

```
# Perform SCT-normalization on each individual sample
seurat_ls <- SplitObject(seurat_clean, split.by = "orig.ident")

for (i in 1:length(seurat_ls)) {
  print(paste0("Processing sample ", i))
  seurat_ls[[i]] <- SCTransform(seurat_ls[[i]], verbose = TRUE,
                                variable.features.n = 3000,
                                vars.to.regress = rna_regress)
}

# Find integration features across samples
integ_features <- SelectIntegrationFeatures(object.list = seurat_ls, nfeatures = 3200) 

# Exclude TRAV/TRAJ and TRBV/TRBJ genes from top features for integration
delIdx <- c(grep("TRBV", integ_features), grep("TRBJ", integ_features), 
            grep("TRAV", integ_features), grep("TRAJ", integ_features))
sort(integ_features[delIdx])
integ_features <- integ_features[-delIdx]
integ_features <- integ_features[1:3000]

# Prepare SCT objects for integration
seurat_ls <- PrepSCTIntegration(object.list = seurat_ls,
				anchor.features = integ_features)

# Find anchors across normalized dataset
integ_anchors <- FindIntegrationAnchors(object.list = seurat_ls, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)

# Integrate by sample
seurat_cca <- IntegrateData(anchorset = integ_anchors,
                            normalization.method = "SCT")

# Update UMAP using integrated assay
DefaultAssay(seurat_cca) <- "integrated"
seurat_cca <- RunPCA(seurat_cca)
seurat_cca <- RunUMAP(seurat_cca, dims = 1:40)
```

See normalization and harmonization reports in [this GitHub repo](https://github.com/hbc/sinclair_scRNA-CITE-seq_mouse_T-cells_aging_hbc04771/tree/main) (single sample) or [this one](https://github.com/hbc/chandraker_scRNASeq_human_PBMC_kidney_transplant_hbc04749/tree/main) (multiple samples with merge/integration) for a complete example.

***
<br>
<br>



### Considerations and Best Practices
<br>
<br>

Whether or not to remove TCR and BCR transcripts will be experiment-specific. Always consult with the client about these decisions. In general, an iterative approach is probably warranted, by initially observing clustering and annotation with the transcripts retained in the object. This can help with annotation, help to potentially identify biologically interesting phenomena like lymphocyte clonal expansion or differences in BCR light chain ratios, and help to identify and annotate rare or novel cell subsets.

Only remove transcripts to observe their influence on clustering or to remove "nuisance" clustering that is not experimentally or biologically interesting.

If you are removing transcripts, always double check the names of the transcripts you are targeting. It is easy to mistakenly include non-antigen receptor transcripts that have similar lettering in the symbols. 

Where best in the workflow to filter antigen receptor features? I'd advocate making a new filtered Seurat Object for comparison. And I would re-do at least the clustering workflow (for example PCA > UMAP > Harmony > Harmony UMAP) using the filtered RNA data, and you may want to renormalize the data as well.

You may want to remove transcripts from both the RNA and SCT assays, which will require two different indexing vectors taken from the rownames of these different assays.

An alternative approach, see above Option 2.

<br>
<br>

**MAIT cells and iNKT cells**

[MAIT cells](https://en.wikipedia.org/wiki/Mucosal-associated_invariant_T_cell)



[iNKT cells](https://en.wikipedia.org/wiki/Natural_killer_T_cell)



Mucosal Associated Invariant T (MAIT) cells and Invariant Natural Killer-like T (iNKT) cells are a subset of T cells, or T cell-like cells that utilize a limited repertoire of TCR gene segments, and the presence of these transcripts could potentially be useful for identifying these rare cell subsets.

In humans, MAIT cells have been documented as using TRAV1-2 (official nomenclature V*a*7.2) and TRAJ12, TRAJ20, and TRAJ33. (Mice use a Trav1-Traj33 TCR*a* chain, paired with Trbv19 or Trbv13 TCR*b* chains.

(Importantly, look for MAIT-specific antibodies in CITE-Seq panels ("TCRVa7.2-prot" for example, in the Biologends TotalSeq B panel. Also note the difference in nomenclature: MAIT-specific _gene_ alpha chain = TRAV**1-2**; MAIT-specific _protein_ alpha chain = TCRVa**7.2**...) MAIT cells may best be best found using integrated RNA-Seq and ADT signals.)


Human iNKT cells have been documented as using TRV24, TRAJ18, and TRBV11, (mice use Trav14 and Traj18).

However, TCR usage in these cells types is an ongoing avenue of research, so it may be better to identify and annotate these cell types using unfiltered datasets and other known markers.

<br>
<br>

**Vertebrate species with poorly annotated transcripts**

Some vertebrate species may not have sufficient annotation to allow identification of TCR and BCR transcripts by Gene Symbols. The IMGT provides correspondence between species [here](https://www.imgt.org/IMGTrepertoire/LocusGenes/#N).

This code, (provided by Will Gammerdinger, thanks Will!!) should help identify these transcripts in GFF3 annotation.

````
This one-liner should give you a list of transcript IDs from a GFF3 that are labeled as V/J/C gene segments:
awk '$3 ~ /[VJC]_gene_segment/' your_organism.gff3 | awk '{print $9}' | awk -F ";" '{print $1}' | sed 's/ID=//' 
````










