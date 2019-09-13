---
title: Single Cell Quality control with saturation
category: Single Cell
---

People often ask how many cells they need to sequence in their next experiment.
Saturation analysis helps to answer that question looking at the current experiment.
Would adding more coverage to the current experiment result in getting more transcripts,
genes, or in just more duplicated reads?

First, use [from bcbio to single cell script](https://github.com/hbc/hbcABC/blob/master/inst/rmarkdown/Rscripts/singlecell/from_bcbio_to_singlecell.R)
to load data from bcbio into SingleCellExperiment object.

Second, use [this Rmd template](https://github.com/hbc/hbcABC/blob/master/inst/rmarkdown/templates/simple_qc_single_cell/skeleton/skeleton.rmd)
to create report.
