---
title: Tools related to RNAseq 
description: This page shows tools that can be applied to RNAseq data.
category: research
subcategory: rnaseq
tags: [literature]
---


# Alternative Splicing

* [IsoformSwitchAnalyzer](https://bioconductor.org/packages/release/bioc/vignettes/IsoformSwitchAnalyzeR/inst/doc/IsoformSwitchAnalyzeR.html) helps to detect alternative splicing. I tried it and an example of a consults is here:https://code.harvard.edu/HSPH/hbc_RNAseq_christiani_RNAediting_on_lung_in_humna_hbc02307. This packages has very nice figures: https://www.dropbox.com/work/HBC%20Team%20Folder%20(1)/consults/david_christiani/RNAseq_christiani_RNAediting_on_lung_in_humna?preview=dtu.html (see at the end of the report).

* [DEXseq]:Following this paper from MLove et al: https://f1000research.com/articles/7-952/v3 I used salmon and DEXseq to call isoform switching. This consult has an example: https://code.harvard.edu/HSPH/hbc_RNAseq_christiani_RNAediting_on_lung_in_humna_hbc02307. I found that normally one isomform changes a lot and another very little, but I found some examples were the switching is more evident.



# Normalization
* [Comparing the normalization methods for
the differential analysis of Illumina high-
throughput RNA-Seq data](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0778-7). Paper that compares different normalization methods on RNASeq data.

# Power
* [Power in pairs: assessing the statistical value of paired samples in tests for differential expression](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6302489/). Paper that looks at effect of paired-design on power in RNA-seq.
