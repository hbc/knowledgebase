---
title: General NGS
description: Quick tips for ngs data.
category: research
subcategory: general_ngs
tags: [ngs]
---

# 3' DGE from LSP demultiplexing example
```
bcl2fastq --adapter-stringency 0.9 --barcode-mismatches 0 --fastq-compression-level 4 --min-log-level INFO --minimum-trimmed-read-length 0 --sample-sheet /n/boslfs/INSTRUMENTS/illumina/180604_NB501677_0276_AHVMT2BGX5/SampleSheet.csv --runfolder-dir /n/boslfs/INSTRUMENTS/illumina/180604_NB501677_0276_AHVMT2BGX5 --output-dir /n/boslfs/ANALYSIS/180604_NB501677_0276_AHVMT2BGX5 --processing-threads 8 --no-lane-splitting --mask-short-adapter-reads 0 --use-bases-mask y*,y*,y*,y*
```

# Basespace

Use Illumina's native GUI client or run [BaseMount](https://basemount.basespace.illumina.com) on Ubuntu.

The [Python downloader](https://support.basespace.illumina.com/knowledgebase/articles/403618-python-run-downloader) is deprecated and no longer supported by Illumina.



# Figuring out what sequencer was used from FASTQ read names

(via stackoverflow)

```
AAXX = Genome Analyzer
BCXX = HiSeq v1.5
ACXX = HiSeq High-Output v3
ANXX = HiSeq High-Output v4
ADXX = HiSeq RR v1
AMXX, BCXX =HiSeq RR v2
ALXX = HiSeqX
BGXX, AGXX = High-Output NextSeq
AFXX = Mid-Output NextSeq
5 letter/number = MiSeq
```
