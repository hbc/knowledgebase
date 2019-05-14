---
title: Tools related to ChIP-seq and related analysis 
description: This page shows tools that can be applied to ChIP-seq analysis and other chromatin biology related techniques.
category: research
subcategory: chipseq 
tags: [tools, literature]
---

## Tools_tested

### Alignment

The most **commonly used aligner for ChIP-seq is Bowtie2**. However, recently we have tried runing bwa which has resulted in higher mapping rates (~ 2%), with an equally similar increase in the number of duplicate mappings identified. Post-filtering this translates to a significantly higher number of mapped reads and results in a much larger number of peaks being called (30% increase). When we compare the peak calls generated from the different aligners, the bwa peak calls are a superset of those called from the Bowtie2 aligments. Whether or not these additional peaks are true positives, is something that is yet to be determined.

- bwa
 - Rory Kirchner, tested in early 2018
    - version# ?
    - information about parameters used?
    - TODO: find a good benchmarking dataset to compare bwa and 
 
### Peak calling

- MACS2
  - this is what we currently use in bcbio. Has both narrow and broad peak functionality and is used for ATAC-seq peak calling.

- SPP
  - Meeta Mistry, tested/used last in 2017
  		- R 3.2.1	
  		- not sure if it works well for broad peaks 
  		- had trouble getting it to work with more recent R versions
  

### QC

- [ChIPQC](http://bioconductor.org/packages/release/bioc/html/ChIPQC.html)
	- currently being used in most consults
	- the report is easy to generate (two lines of code)
	- Run locally with O2 mount or Run on O2??		- locally can be computationally intensive since processing BAM files
		- getting all dependency packages installed on the cluster is a bit of an issue
		-  Solution: Can we have it run as part of bcbio?

	
- [phantompeakqualtools](https://github.com/kundajelab/phantompeakqualtools) 
  - Meeta Mistry, last used in 2018
  - not as relevant for braod peak calling
  - no longer used (or presented in teaching materials) since some of these metrics are covered in ChIPQC

	

### Handling replicates

- bedtools
  - standard, quick way of checking overlapping peaks between replicates
  - min 1bp overlap
  - no statistical significance

- IDR (Irreproducibility Discovery Rate) 
  - a rank-based method of evaulating concordance between peak calls. Takes the list of overlapping peaks (1bp overlap) and statistically evaluates when the concordance drops off. IDR values are assigned to each peak, which can be interepreted similar to an FDR
  - part of the ENCODE guidelines but not sure which version of the tool is best to use
    - original pipeline is [deprecated](https://sites.google.com/site/anshulkundaje/projects/idr)
    - links to a new pipeline [which is also deprecated](https://github.com/kundajelab/chipseq_pipeline)
    - the latest pipeline, which is actually a processing pipeline (not just IDR) is [apparently actively in use](https://github.com/ENCODE-DCC/chip-seq-pipeline2)
    - we have been using [IDR from @nboley](https://github.com/nboley/idr) for consults and teaching (it is a module on O2) but MM recently [submitted an issue](https://github.com/nboley/idr/issues/43) but have not heard back. Not sure if this is being actively maintined or not.
    - there is also an [IDR package in CRAN](https://cran.r-project.org/web/packages/idr/index.html)

### Visualization

- deepTools

### Differential Enrichment

### Functional analysis and Annotation 

- ChIPseeker 
	- R Bioconductor package for peak annotation and visualization. Currently, in use for most 2018 consults onwards. 
	- Nearest gene annotation, uses the TxDb databases.
	- Visualization is based on peaks, not read density - so not very accurate.
	- Target gene lists can be used directly as input to clusterProfiler for functional analysis

- HOMER 
	- Meeta Mistry, used this for two big ChIP-seq consults (Harwell and Flanagan). Last used in 2017.
	- Good for peak annotation and can also do some level of visualization (generated the underlying data which can be loaded into R for plotting). 
	- Internally uses RefSeq by default, but can provide a GTF file for custom annotation. 
	- No functional analysis, but is useful for motif analysis.



## Tools_novel


- [csaw](https://bioconductor.org/packages/release/bioc/html/csaw.html)
  - Detection of differentially bound regions in ChIP-seq data with sliding windows, with methods for normalization and proper FDR control.
 

- [SICER's reimplementation: epic2](https://github.com/biocore-ntnu/epic2)
  - Link to the [paper](https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/btz232/5421513?redirectedFrom=fulltext)
  - It might be worth taking some time to investigate peak callers designed specifically for broad marks. We default to MACS2 `--broad` but depedning on the histone mark have had trouble finding peaks.

- [SUPERmerge](https://www.biorxiv.org/content/10.1101/121897v1)
 - broad peak caller; especially useful for low sample sizes

- [haystack_bio](https://github.com/pinellolab/haystack_bio)
 - An analysis pipeline from the Pinello lab. It can be used with histone modifications and chromatin accessibility data generated by ChIP-seq, DNase-Seq, and ATAC-seq assays across multiple cell-types. In addition, it is also possible to integrate gene expression data generated by RNA-seq for example.

- [TF and Histone ChIP-seq processing pipeline](https://github.com/ENCODE-DCC/chip-seq-pipeline2)
  - from the Kundaje lab


