---
title: Stranded RNA-seq libraries.
description: Explains strandedness and where to find info in bcbio.
category: research
subcategory: rnaseq 
---

Bulk RNA-seq libraries retaining strand information (stranded) are useful to quantify expression with higher accuracy for opposite 
strand transcripts which overlap or have overlapping UTRs.
https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-015-1876-7. 

Bcbio RNA-seq pipeline has a 'strandedness' parameter: [unstranded|firststrand|secondstrand]  
https://bcbio-nextgen.readthedocs.io/en/latest/contents/pipelines.html?highlight=strand#configuration.

The terminology was inherited from Tophat, see the detailed description in the Salmon doc. 
https://salmon.readthedocs.io/en/latest/library_type.html
Note, that firstrand = ISR.

If the strandedness is unknown, run a small subset of reads with 'unstranded' in bcbio and check out what Salmon reports in 
`bcbio_project/final/sample/salmon/lib_format_counts.json`:
```
{
    "read_files": [
        "/dev/fd/63",
        "/dev/fd/62"
    ],
    "expected_format": "IU",
    "compatible_fragment_ratio": 1.0,
    "num_compatible_fragments": 721856,
    "num_assigned_fragments": 721856,
    "num_frags_with_concordant_consistent_mappings": 692049,
    "num_frags_with_inconsistent_or_orphan_mappings": 47441,
    "strand_mapping_bias": 0.9477291347866986,
    "MSF": 0,
    "OSF": 0,
    "ISF": 36174,
    "MSR": 0,
    "OSR": 0,
    "ISR": 655875,
    "SF": 37676,
    "SR": 9765,
    "MU": 0,
    "OU": 0,
    "IU": 0,
    "U": 0
}
```
Here the majority of reads are ISR.

Another way to check strand bias is  
`bcbio_project/final/sample/qc/qualimap_rnaseq/rnaseq_qc_results.txt`.  
It has `SSP estimation (fwd/rev) = 0.04 / 0.96` meaning strand bias (ISR, firststrand).

Yet another way to confirm strand bias is seqc.  
http://rseqc.sourceforge.net/#infer-experiment-py.  
It uses a small subset of the input bam file: 
`infer_experiment.py -r /bcbio/genomes/Mmusculus/mm10/rnaseq/ref-transcripts.bed -i test.bam`

```
This is PairEnd Data
Fraction of reads failed to determine: 0.1461
Fraction of reads explained by "1++,1--,2+-,2-+": 0.0177
Fraction of reads explained by "1+-,1-+,2++,2--": 0.8362
```
