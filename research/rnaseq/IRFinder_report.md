---
title: Intron retention analysis
description: This code helps with the analysis of intron retention from RNAseq data.
category: research
subcategory: rnaseq
tags: [intron_retention, differential_expression]
---


## Methods for IRFinder

The intron retention analysis was performed using IRFinder version 1.2.3 based on the gene annotations for the hg19 reference genome.

- All potential introns were analyzed for intron retention ratios in each condition, and introns were defined as any region between two exon features in any transcript. The exon coordinates were derived from reference gene annotation files.

- Intronic regions with poor mappability based on pre-determined mapping of synthetic reads to the genome were excluded from the measurable intron area.

- Intron retention was assessed for directional RNA-Seq data (reverse strand).

- The Audic and Claverie test for small numbers of replicates (Max = 3 replicates) was applied to the quantified intron retention events to determine differential expression between Control and NRDE2-KD conditions as described in [this post](http://mimirna.centenary.org.au/irfinder/example1.html).

For more detailed descriptions, see the [IRFinder wiki](https://github.com/williamritchie/IRFinder/wiki).

```r
library(tidyverse)
library(annotables)

# Small Amounts of Replicates via Audic and Claverie Test
irfinder_AC_test <- read.table("/path/to/IRFinder/results.tab", header = T)

parsed_rownames <- str_split(irfinder_AC_test$Intron.GeneName.GeneID, "/", simplify = TRUE)

irfinder_AC_test <- data.frame(irfinder_AC_test, parsed_rownames)

irfinder_AC_test <- irfinder_AC_test %>% rename(gene=X1,
           ensembl_id=X2,
           splicing=X3)
grch37_description <- grch37[, c("ensgene", "biotype", "description")]

grch37_description <- grch37_description[which(!(duplicated(grch37_description$ensgene))), ]

irfinder_ACtest_merged <- merge(irfinder_AC_test, grch37_description, by.x="ensembl_id", by.y="ensgene")

irfinder_ACtest_merged <- irfinder_ACtest_merged[, -5]

# Need to adjust for multiple test correction
irfinder_ACtest_merged$padj <- p.adjust(irfinder_ACtest_merged$p.diff, "BH")

# Order by padj
irfinder_ACtest_merged <- irfinder_ACtest_merged[order(irfinder_ACtest_merged$padj), ]

# Write to file all results
write.csv(irfinder_ACtest_merged, "results/irfinder_ACtest_all_results_padj.csv")

# Determine significant results with padj < 0.05
sig_irfinder_ACtest <- irfinder_ACtest_merged[which(irfinder_ACtest_merged$padj < 0.05), ]

sig_irfinder_ACtest <- sig_irfinder_ACtest[order(sig_irfinder_ACtest$p.diff),]
```

## Results

There were ### significantly retained introns. The results output for each intron includes the following information, also described in the [IRFinder wiki](https://github.com/williamritchie/IRFinder/wiki) and an [analysis example](http://mimirna.centenary.org.au/irfinder/example1.html).

- **ensembl_id:** Ensembl ID
- **Chr**: chromosome
- **Start:** start coordinates of intron
- **End:** end coordinates of intron
- **Direction:** strand
- **ExcludedBases:** number of bases within the intronic region that have been excluded from the calculation of intronic coverage because of overlapping features or mapping issues. The 5bp flanking any exon are also excluded, hence all introns exclude at least 10 bases.
- **p.diff:** p-value output from the Audic and Claverie test regarding intron retention differences between conditions. Indicates to what extent intron-retention has significantly changed between Condition A (control) and Condition B (Nrde2-KD)
- **p.increased:** p-value for whether higher intron retention in Nrde2-KD relative to control
- **p.decreased:** p-value for whether lower intron retention in Nrde2-KD relative to control
- **A.IRratio:** intron retention ratio for control. Calculated by: IntronDepth / (max(number of reads that map the 3' flanking exon surrounding the intron and to another exon within the same gene, number of reads that map the 5' flanking exon surrounding the intron and to another exon within the same gene) + IntronDepth)
- **A.IRok:** Warnings about potential biases to the IR calculation. Low coverage, non-uniform coverage of the intron, etc. "NonUniformCover" indicates that the multiple places IRFinder measures the depth of the intron are not consistent. A visual check of the RNA-Seq trace is advised to rule out alternate-TSS or similar. In some cases this warning may simply be triggered by the uneven cover often seen in short-read RNA-Seq experiments.
- **A.IntronCover:** Ratio of bases with mapped reads for control
- **A.IntronDepth:** number of reads that map over a given bp for control. IntronDepth is the median depth of the intronic region without the excluded regions. It is used to calculate the IRratio. Excluded regions comprise ExclBases and bases with the top and bottom 30% of intronic depth.
- **A.SplicesMax:** number of reads that map to any exon for control? (no documentation)
- **A.SplicesExact:** number of reads that map across the 3' and 5' flanking exons for control
- **B.IRratio:** intron retention ratio for Nrde-KD.
- **B.IRok:** Warnings about potential biases to the IR calculation for Nrde-KD.
- **B.IntronCover:** Ratio of bases with mapped for Nrde-KD.
- **B.IntronDepth:** number of reads that map over a given bp for Nrde-KD.
- **B.SplicesMax:** number of reads that map to any exon for Nrde-KD? (no documentation).
- **B.SplicesExact:** number of reads that map across the 3' and 5' flanking exons for Nrde-KD.
- **replicates:** whether or not there were replicates (no documentation)
- **A1.IRratio:** IR ratio for first replicate for control condition
- **A2.IRratio:** IR ratio for second replicate for control condition
- **A3.IRratio:** IR ratio for third replicate for control condition
- **B1.IRratio:** IR ratio for first replicate for Nrde2-KD condition
- **B2.IRratio:** IR ratio for second replicate for Nrde2-KD condition
- **B3.IRratio:** IR ratio for third replicate for Nrde2-KD condition
- **gene:** associated gene name
- **splicing:** whether there are known transcripts that include this intron. `known-exon` included if there are known transcripts.
- **biotype:** type of transcript (protein_coding, mtRNA, etc.)
- **description:** short description of gene
- **padj:** p-value corrected for multiple testing using the Benjamini-Hochberg/FDR method

Similar to the criteria used in the [IRFinder paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1184-4), the significant intron retention events were defined as:

- `p.adj` < 0.05
- `IRratio` for either control or treatment > 0.1 (at least 10% of transcripts exhibit retention of the intron for at least one of the conditions)
- Retained introns exhibit a coverage (`IntronDepth`) of at least three reads across the entire intron after excluding non-measurable intronic regions

All results can be accessed using the links below the table. The top 20 most significant retained introns are displayed below:

```r
# Only returning those introns represented in more than 10% of transcripts in A or B
sig_irfinder_filtered <- sig_irfinder_ACtest[which(sig_irfinder_ACtest$A.IRratio > 0.1 | sig_irfinder_ACtest$B.IRratio > 0.1), ]

# Only returning those introns with a coverage of more than three reads across the entire intron after excluding non-measurable intronic regions
sig_irfinder_filtered2 <- sig_irfinder_filtered[sig_irfinder_filtered$A.IntronDepth > 3 | sig_irfinder_filtered$B.IntronDepth > 3, ]

sig_irfinder_filtered2 <- sig_irfinder_filtered2[order(sig_irfinder_filtered2$padj), ]
sig_irfinder_filtered2$p.diff <- formatC(sig_irfinder_filtered2$p.diff, format = "e", digits = 2)
sig_irfinder_filtered2$p.increased <- formatC(sig_irfinder_filtered2$p.increased, format = "e", digits = 2)
sig_irfinder_filtered2$p.decreased <- formatC(sig_irfinder_filtered2$p.decreased, format = "e", digits = 2)
knitr::kable(sig_irfinder_filtered2[1:20, ])
write.csv(sig_irfinder_filtered2, "results/significant_irfinder_ACtest_results_padj.csv")
```

[Download all results]

[Download significant results]

The significantly retained introns were explored for several genes of interest.

```r
interesting_genes <- sig_irfinder_filtered2[sig_irfinder_filtered2$gene %in% c(), ]
knitr::kable(interesting_genes)

```

```r
sessionInfo()
```
