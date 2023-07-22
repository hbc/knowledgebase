## Musings on T Cell Receptor (TCR) and B Cell Receptor (BCR) V,D,J and C transcripts, and their influence on clustering in single-cell RNA-Seq experiments
<br>
<br>

In vertebrates, the lymphocyte antigen receptors, T Cell Receptor (TCR) and B Cell Receptor (BCR) are encoded by V ("variable"), D ("diversity"), J("joining"), and C ("constant") genes (or gene "segments") located in TCR and BCR loci. Numerous different versions of each segment type are found in germline DNA. (For example, in humans there are > 80 different TCR*a* V segments.) During lymphocyte development, through a process called "somatic recombination", or "v(d)j recombination" (in contrast to "genetic recombination" during meiosis), individual v,d,j and c segments are selected to encode a unique TCR or BCR sequence for each cell. The unselected segments are excised from the DNA. The nearly random selection and combination of individual v,d,j and c segments, promotes a vast (est. > 10<sup>15</sup> ) combinatorial diversity of potential unique antigen binding domains encoded on individual lymphocyte antigen receptors. Antigen receptor diversity generated through somatic recombination, and immune memory promoted through lymphocyte clonal expansion are the hallmarks of vertebrate adaptive immunity.


[https://en.wikipedia.org/wiki/V(D)J_recombination](https://en.wikipedia.org/wiki/V\(D\)J_recombination)

[v(d)j recombination](https://www.ncbi.nlm.nih.gov/books/NBK27140/)
<br>
<br>

**Antigen receptor transcripts and clustering**


The presence of RNA transcript counts arising from antigen receptor sequences can promote clustering in sc-RNA-Seq analysis. This can be biologically interesting, or not, depending on the experiment.

For example, there are two subtypes of TCRs, TCR*a*/*b*, or TCR*g*/*d*.  These two different T cell subsets  have distinct biological activities, so segregation of these two subsets is generally useful, and the differential expression of TRA, TRB, TRG and TRD transcripts are useful to promote subset segregation and annotation.

There are also two general subtypes of BCRs, defined by the usage of either lambda ( *k* ) or kappa ( *l* ) light chains. The biological activities, however, of these two B cell subtypes are less distinct, so clustering induced by differential *l* or *k* gene segment usage may not be interesting to a specific experiment, and in fact, can often promote significant subclustering which can interfere with identification and annotation of more biologically or experimentally interesting clustering and annotation. (Although, there may be some experimental situations where quantification of kappa and lambda B cells is of interest.)
<br>
<br>


**Removing antigen receptor transcripts from Seurat objects**


If the analyst and the client are interested in observing the effects on clustering promoted by antigen receptor genes, these are readily removed by subsetting the object.

This code displays the names of the TCR*a* chain V and J transcripts found in the Seurat object. (TCR*a* and *d* chains don't utilize D segments)

````
DefaultAssay(SeuratObject)<-"RNA"

rownames(SeuratObject)[(grep(pattern = '^TRAV|^TRAJ' , rownames(SeuratObject)))]
````
		
which will return something like this:
````
[1] "TRBV1"       "TRBV2"       "TRBV3-1"     "TRBV4-1"     "TRBV5-1"    
[6] "TRBV6-1"     "TRBV7-1"     "TRBV4-2"     "TRBV6-2"     "TRBV7-2"    
[11] "TRBV8-1"     "TRBV5-2"     "TRBV6-4"     "TRBV7-3"     "TRBV8-2"    
[16] "TRBV5-3"     "TRBV9"       "TRBV10-1"    "TRBV11-1"    "TRBV12-1"   
[21] "TRBV10-2"    "TRBV11-2"    "TRBV12-2"    "TRBV6-5"     "TRBV7-4"    
[26] "TRBV5-4"     "TRBV6-6"     "TRBV7-5"     "TRBV5-5"     "TRBV6-7"    
[31] "TRBV7-6"     "TRBV5-6"     "TRBV6-8"     "TRBV7-7"     "TRBV5-7"    
[36] "TRBV7-9"     "TRBV13"      "TRBV10-3"    "TRBV11-3"    "TRBV12-3"   
[41] "TRBV12-4"    "TRBV12-5"    "TRBV14"      "TRBV15"      "TRBV16"     
[46] "TRBV17"      "TRBV18"      "TRBV19"      "TRBV20-1"    "TRBV21-1"   
[51] "TRBV22-1"    "TRBV23-1"    "TRBV24-1"    "TRBV25-1"    "TRBVA"      
[56] "TRBV26"      "TRBVB"       "TRBV27"      "TRBV28"      "TRBV29-1"   
[61] "TRBJ1-1"     "TRBJ1-2"     "TRBJ1-3"     "TRBJ1-4"     "TRBJ1-5"    
[66] "TRBJ1-6"     "TRBJ2-1"     "TRBJ2-2"     "TRBJ2-2P"    "TRBJ2-3"    
[71] "TRBJ2-4"     "TRBJ2-5"     "TRBJ2-6"     "TRBJ2-7"     "TRBV30"     
[76] "TRBV20OR9-2" "TRBV21OR9-2" "TRBV22OR9-2" "TRBV23OR9-2" "TRBV24OR9-2"
[81] "TRBV25OR9-2" "TRBV26OR9-2" "TRBV29OR9-2"
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


**Considerations and Best Practices**
<br>
<br>

Whether or not to remove TCR and BCR transcripts will be experiment-specific. Always consult with the client about these decisions. In general, an iterative approach is probably warranted, by initially observing clustering and annotation with the transcripts retained in the object. This can help with annotation, help to potentially identify biologically interesting phenomena like lymphocyte clonal expansion or differences in BCR light chain ratios, and help to identify and annotate rare or novel cell subsets.

Only remove transcripts to observe their influence on clustering or to remove "nuisance" clustering that is not experimentally or biologically interesting.

If you are removing transcripts, always double check the names of the transcripts you are targeting. It is easy to mistakenly include non-antigen receptor transcripts that have similar lettering in the symbols. 

Where best in the workflow to filter antigen receptor features? I'd advocate making a new filtered Seurat Object for comparison. And I would re-do at least the clustering workflow (for example PCA > UMAP > Harmony > Harmony UMAP) using the filtered RNA data, and you may want to renormalize the data as well.

You may want to remove transcripts from both the RNA and SCT assays, which will require two different indexing vectors taken from the rownames of these different assays.

An alternative approach is just to filter receptor transcripts from the scale data slot in the SCT assay.

<br>
<br>

**MAIT cells and iNKT cells**

[iNKT cells](https://en.wikipedia.org/wiki/Natural_killer_T_cell)
[MAIT cells](https://en.wikipedia.org/wiki/Mucosal-associated_invariant_T_cell)


Mucosal Associated Invariant T (MAIT) cells and Invariant Natural Killer-like T (iNKT) cells are a subset of T cells, or T cell-like cells that utilize a limited repertoire of TCR gene segments, and the presence of these transcripts could potentially be useful for identifying these rare cell subsets.

In humans, MAIT cells have been documented as using TRAV1-2 and TRAJ12, TRAJ20, and TRAJ33. (Mice use a Trav1-Traj33 TCR*a* chain, paired with Trbv19 or Trbv13 TCR*b* chains.

(Importantly, look for MAIT-specific antibodies in CITE-Seq panels ("TCRVa7.2-prot" for example, in the Biologends TotalSeq B panel.) MAIT cells may best be best found using integrated RNA-Seq and ADT signals.)


Human iNKT cells have been documented as using TRV24, TRAJ18, and TRBV11, (mice use Trav14 and Traj18).

However, TCR usage in these cells types is an ongoing avenue of research, so it may be better to identify and annotate these cell types using unfiltered datasets.

<br>
<br>

**Vertebrate species with poorly annotated transcripts**

Some vertebrate species may not have sufficient annotation to allow identification of TCR and BCR transcripts by Gene Symbols.

This code, (provided by Will Gammerdinger, thanks Will!!) should help identify these transcripts in GFF3 annotation.

````
This one-liner should give you a list of transcript IDs from a GFF3 that are labeled as V/J/C gene segments:
awk '$3 ~ /[VJC]_gene_segment/' your_organism.gff3 | awk '{print $9}' | awk -F ";" '{print $1}' | sed 's/ID=//' 
````










