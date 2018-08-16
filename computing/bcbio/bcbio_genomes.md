---
title: Genomes in bcbio
description: This describes redundant genomes in bcbio and their differences
category: computing
subcategory: bcbio
tags: [bcbio]
---

**Drosophila melanogaster**

*DGP6* - built from Flybase
  - has a different format for annotations for non-coding genes in the gtf
  - only protein coding genes will make it into Salmon and downstream
  
*DGP6.92* - built from Ensembl info
  - will have all non-coding RNAs in Salmon and downstream results
  - shows lower gene detection rates than Flybase
