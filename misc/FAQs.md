# Snippets for Frequently Asked Questions from clients, when they go over the reports.
```
Feel free to add the FAQs you have received.
```

## General
### Functional Analysis
1. What is geneRatio and bgRatio in overerpresentation analysis?

- The geneRatio is the {# of annotated genes assigned to term from input}/{# of input genes annotated}

- The bgRatio is the {# of annotated genes assigned to term from background}/{# of background genes annotated}

Please note that the denominator may be different between MF, BP, and CC as there are different number of genes annotated for those categories.

- The input is a list of candidate genes (i.e. list of significant DEGs) while the background is a list of all the genes in the study.

2. How is the p-value calculated?

Simplistic link: https://www.pathwaycommons.org/guide/primers/statistics/fishers_exact_test/

A bit more mathematical link: http://www.nonlinear.com/progenesis/qi/v2.0/faq/should-i-use-enrichment-or-over-representation-analysis-for-pathways-data.aspx

Good video link: https://www.coursera.org/lecture/bd2k-lincs/enrichment-analysis-part-1-xLgN5

### Multiple testing correction
What is q-value and why do we need this?

Here is a pretty good slide about the need for multiple testing and how FDR is calculated : 
https://www.gs.washington.edu/academics/courses/akey/56008/lecture/lecture10.pdf

### t-test and Wilcoxon {credits to [Preeti](https://github.com/orgs/hbc/people/preetida)}
Blog by Jonathan Bartlett, super helpful for many stat related questions.

https://thestatsgeek.com/2014/04/12/is-the-wilcoxon-mann-whitney-test-a-good-non-parametric-alternative-to-the-t-test/

## ChIP-seq and ATAC-seq

## BULK RNA-seq

## scRNA-seq

