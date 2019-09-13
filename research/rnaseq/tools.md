---
title: Tools related to bulk RNA-seq analysis 
description: This page shows tools that can be applied to bulk RNA-seq analysis.
category: research
subcategory: rnaseq 
tags: [tools, literature]
---

## Tools_tested

> - Exampletool_name (as link to tool)
>     - tester initials, date of testing (yyyy_mm)
>     - version#
>     - brief description about what the tool does
>     - information about how the tool performed
>     - what requirements are needed (e.g. R-3.5.1, etc.)
>     - whether we have teaching material for this tool with links
>     - whether incorporated in bcbio or downstream analysis template/report (if so, then include date of incorporation)
>     - if used in a template/report, then provide a link to the template/report
    
 - [IsoformSwitchAnalyzer](https://bioconductor.org/packages/release/bioc/vignettes/IsoformSwitchAnalyzeR/inst/doc/IsoformSwitchAnalyzeR.html)
    - LP/VB?, 2019/02?
    - version#?
    - helps to detect alternative splicing
    - output very nice figures
    - what requirements are needed (e.g. R-3.5.1, etc.)?
    - no tutorials available
    - not incorporated into bcbio
    -  I tried it and an example of a consults is here:https://code.harvard.edu/HSPH/hbc_RNAseq_christiani_RNAediting_on_lung_in_humna_hbc02307. This packages has very nice figures: https://www.dropbox.com/work/HBC%20Team%20Folder%20(1)/consults/david_christiani/RNAseq_christiani_RNAediting_on_lung_in_humna?preview=dtu.html (see at the end of the report).
    
- [DEXseq](https://bioconductor.riken.jp/packages/3.0/bioc/html/DEXSeq.html)
    - LP/VB/RK?, date?
    - version#?
    - used to call isoform switching
    - not recommended - use DTU tool instead
    - what requirements are needed (e.g. R-3.5.1, etc.)?
    - no tutorials available
    - yes, DEXseq is incorporated in bcbio
    - Following this paper from MLove et al: https://f1000research.com/articles/7-952/v3 I used salmon and DEXseq to call isoform switching. This consult has an example: https://code.harvard.edu/HSPH/hbc_RNAseq_christiani_RNAediting_on_lung_in_humna_hbc02307. I found that normally one isomform changes a lot and another very little, but I found some examples were the switching is more evident.



## Tools_novel

> - Exampletool_name:
> - brief description of tool (one-line)
>   - link to tool
>   - initials of person planning to test tool
