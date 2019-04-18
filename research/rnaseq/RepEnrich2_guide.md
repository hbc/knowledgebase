---
title: How to run repeat enrichment analysis
description: This guide shows how to run RepEnrich2
category: research
subcategory: rnaseq
tags: []
---

RepEnrich2 tries to look at something that standard RNA-seq pipelines miss, the
enrichment of repeats in NGS data. It is extremely slow and is a pain to get
going. Below is a guide getting it working and has some links to a fork of
RepEnrich2 I made that makes it more friendly to use.

I have not actually validated the RepEnrich2 output, so caveat emptor.

# Preparing RepEnrich2

## Create isolated conda environment

```bash
conda create -c bioconda -n repenrich2 python=2.7 biopython bedtools samtools bowtie2 bcbio-nextgen
```

## Download my fork of RepEnrich2 
This has quality of life fixes such as memoization of outputs so if it fails you don't
have to redo steps.

```bash
git clone git@github.com:nerettilab/RepEnrich2.git
```

## Download a pre-created index 
You can make your own, for example I made
[hg38](https://www.dropbox.com/s/lefkk38q6bbj76b/Repenrich2_setup_hg38.tar.gz?dl=1)
and the RepEnrich2 folks have mm9 and hg19
[here](https://drive.google.com/drive/folders/0B8_2gE04f4QWNmdpWlhaWEYwaHM). But the RepeatMasker
file it uses needs to be cleaned first and I'm not sure how they cleaned it. They had a hg38 one cleaned
already from RepEnrich so I just used that.

## Download bcbio_RepEnrich2
Download [bcbio_RepEnrich2](https://github.com/roryk/bcbio_RepEnrich2). This will need modification if you
want to use it, but it is simple, I just didn't bother as I don't anticipate us running this again.

# Running RepEnrich2
`bcbio_RepEnrich2` is all you need to run it, the help should give you enough information to go on.
annotation here is the file from RepeatMasker that was used to generate the RepEnrich setup. The
bowtie index is a bowtie2 index of the genome you aligned to. Running RepEnrich2 takes FOREVER, so
be sure to run it on the long queue.

Example command:

```bash
python bcbio_RepEnrich2.py --threads 16 ../human-dsrna/config/human-dsrna.yaml /n/app/bcbio/biodata/genomes/Hsapiens/hg38/bowtie2/hg38 metadata/hg38_repeatmasker_clean.txt metadata/RepEnrich2_setup_hg38/
```

# RepEnrich2 outputs

You will get three files for each sample, for example:

```
P1722_class_fraction_counts.txt
P1722_family_fraction_counts.txt
P1722_fraction_counts.txt
```

The `class` and `family` files are the counts in the `samplename_fraction_counts.txt` file aggregated by family or 
class. Those could be used as aggregate analyses, but the `fraciton_counts` looks at the different repeat
types individually, so is more what folks are probably looking for.
