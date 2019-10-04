---
title: bcbio with multiple lanes
description: This code helps with configuring bcbio to run on samples that have multiple files/lanes
category: computing
subcategory: bcbio
tags: [bcbio]
---

## Configuring bcbio: multiple files per sample

Sequencing facilities often run each sample across multiple lanes, producing technical replicates. Running bcbio normally treats these replicates independent samples when calculating counts. To avoid this, and merge the files correctly prior to counting, we can use the `bcbio_prepare_samples.py` script.

### Prepare a metadata file

Before running this, you need to prepare a metadata file much like the one used in normal bcbio configuration, except containing a row for every file (even paired read files). The samplename column must contain the **full file name**, including **lane number**, **read pair number**, and **extension**. The description column must only be unique at the **sample** level, not the file level, ie: all those files (pairs and lanes) which stemmed from the same sample must have the same description. Here is a **paired end** example:

```
samplename,description,condition
sample1_L001_R1.fastq,s1,wildtype
sample1_L001_R2.fastq,s1,wildtype
sample1_L002_R1.fastq,s1,wildtype
sample1_L002_R2.fastq,s1,wildtype
sample2_L001_R1.fastq,s2,knockout
sample2_L001_R2.fastq,s2,knockout
sample2_L002_R1.fastq,s2,knockout
sample2_L002_R2.fastq,s2,knockout
```

If reads are **single end** it will look more like this:

```
samplename,description,condition
sample1_L001.fastq,s1,wildtype
sample1_L002.fastq,s1,wildtype
sample2_L001.fastq,s2,knockout
sample2_L002.fastq,s2,knockout
```

### Prepare merged samples

Run `bcbio_prepare_samples.py` with flags:   
--out as the output folder where the merged files will land   
--csv as the file-specific metadata prepared above.

```
bcbio_prepare_samples.py --out merged --csv project1.csv
```

The script merges the fastq/bam files to one per sample, with the filename as the corresponding description (eg: `s1.fastq`, `s2.fastq`). If reads are paired, you will get `s1_R1.fastq`, `s1_R2.fastq`, etc. The merged files are in your --out folder, and a new csv `*-merged.csv` is in the same folder as the original csv.

### Prepare bcbio

Now run the `bcbio_nextgen.py` script as normal, but passing it the newly generated csv (`*-merged.csv`), and the new merged fastq files. The yaml remains as normal.

```
bcbio_nextgen.py -w template project1/config/project1-template.yaml project1-merged.csv merged/*fastq
```
   
   
A message from **The Department of Redundancy Department**: This information is more or less repeated [here](https://bcbio-nextgen.readthedocs.io/en/latest/contents/configuration.html#multiple-files-per-sample)











