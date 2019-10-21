---
title: Genomes in bcbio
description: This describes redundant genomes in bcbio and their differences
category: computing
subcategory: bcbio
tags: [bcbio]
---

**Caenorhabditis elegans**

*WBcel235_WS272* - built from wormbase
  - Assembly WBcel235, relase WS272. Project PRJNA13578 (N2 strain)
  - Files:
      - Genomic sequence: ftp://ftp.wormbase.org/pub/wormbase/releases/WS272/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS272.genomic.fa.gz

      - Annotation: ftp://ftp.wormbase.org/pub/wormbase/releases/WS272/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS272.canonical_geneset.gtf.gz

**Drosophila melanogaster**

*DGP6* - built from Flybase
  - has a different format for annotations for non-coding genes in the gtf
  - only protein coding genes will make it into Salmon and downstream
  
*DGP6.92* - built from Ensembl info
  - will have all non-coding RNAs in Salmon and downstream results
  - shows lower gene detection rates than Flybase
 
 **Updating supported transcriptomes**
1. clone cloudbiolinux
2. update transcriptome
```bash
bcbio_python cloudbiolinux/utils/prepare_tx_gff.py --cores 8 --gtf Macaca_mulatta.Mmul_8.0.1.95.chr.gtf.gz --fasta /n/app/bcbio/biodata/genomes/Mmulatta/mmul8noscaffold/seq/mmul8noscaffold.fa Mmulatta mmul8noscaffold
```
3. upload the xz file to the bucket
```bash
aws s3 cp hg19-rnaseq-2019-02-28_75.tar.xz s3://biodata/annotation/ --grants read=uri=http://acs.amazonaws.com/groups/global/AllUsers full=emailaddress=chapmanb@50mail.com
```
4. edit cloudbiolinux ggd transcripts.yaml recipe to point to the new file uploaded on the bucket
5. edit the cloudbiolinux ggd gtf.yaml to show where you got the GTF from and what you did to it
6. test before pushing
```bash
   mkdir tmpbcbio-install
   ln -s `pwd`/cloudbiolinux tmpbcbio-install/cloudbiolinux
   log into bcbio user: sudo -su bcbio /bin/bash
   bcbio_nextgen.py upgrade --data
```
7. push changes back to cloudbiolinux

