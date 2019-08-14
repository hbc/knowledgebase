---
title: Stranded RNA-seq libraries.
description: Explains strandedness and where to find info in bcbio.
category: research
subcategory: rnaseq 
---

Bulk RNA-seq libraries retaining strand information (stranded) are useful to quantify expression with higher accuracy for opposite 
strand transcripts which overlap or have overlapping UTRs (https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-015-1876-7). 

Bcbio RNA-seq pipeline has a 'strandedness' parameter: [unstranded|firststrand|secondstrand]
(https://bcbio-nextgen.readthedocs.io/en/latest/contents/pipelines.html?highlight=strand#configuration).

The terminology was inherited from Tophat, see the detailed description in the Salmon doc (https://salmon.readthedocs.io/en/latest/library_type.html)
Note that firstrand = ISR.

If the strandedness is unknown, run a small portion of reads with unstranded and check out `bcbio_project/final/sample/qc/qualimap_rnaseq/rnaseq_qc_results.txt`.
It has `SSP estimation (fwd/rev) = 0.04 / 0.96` meaning strand bias (ISR, firststrand).

Another way to confirm strand bias is seqc (http://rseqc.sourceforge.net/#infer-experiment-py). It uses a small subset of the input bam file:
`infer_experiment.py -r /bcbio/genomes/Mmusculus/mm10/rnaseq/ref-transcripts.bed -i test.bam`

```
This is PairEnd Data
Fraction of reads failed to determine: 0.1461
Fraction of reads explained by "1++,1--,2+-,2-+": 0.0177
Fraction of reads explained by "1+-,1-+,2++,2--": 0.8362
```
