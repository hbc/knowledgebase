---
title: Improvements for the analysis
description: List of things to try.
category: research
subcategory: orphans
tags: [hbc]
---


1. Try out Alevin from Salmon for a more principled single-cell quantification (https://www.biorxiv.org/content/early/2018/06/01/335000)
2. Add retained intron analysis with IRFinder to bcbio-nextgen
3. See if adding support for grolar to convert pizzly output to something more parseable makes sense. It's an R script and hasn't really been worked on so might not be useable: https://github.com/MattBashton/grolar
4. Add automatic loading/QC of bcbioSingleCell data from bcbio
5. Convert bcbio-nextgen singlecell matrices to HDF5 format in bcbio
6. Swap bcbioSingleCell to read the already-combined matrices for speed purposes
7. Add bcbioRNASeq template to do DTU usage using DRIMseq (https://f1000research.com/articles/7-952/v1)
8. Update installed genomes to use newest Ensembl build for RNA-seq for bcbio-supported genomes.
