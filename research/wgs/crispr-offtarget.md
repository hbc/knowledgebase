---
title: How to find off-target CRISPR edits
description: How to find off-target CRISPR edits
category: research
subcategory: wgs
tags: [annotation]
---

# Overview
This guide is how to call offtarget edits in a CRISPR edited genome. This is
pretty easy to and only takes a few steps. First, we need to figure out what is
different between the CRISPR edited samples and the (hopefully they gave you
these) control samples. Then we need to find a set of predicted off-target
CRISPR sites. Finally, once we know what is different, we need to overlap the
differences with predicted off-target sites, allowing some mismatches.  Then we
can report the overall differences and differences that could be due to
offtarget edits.

# Call CRISPR-edited specific variants
You want to call edits that are in the CRISPRed sample but not the unedited
sample. You can do that by plugging into the tumor-normal calling part of bcbio
and pretending the CRISPR-edited sample is a tumor sample and the non-edited
sample is a normal sample.

To get tumor-normal calling to work you need to use a variant caller that
can handle that, I recommend mutect2.

To tell bcbio that a pair of samples is a tumor-normal pair you need to

1. Put the tumor and normal sample in the same **batch** by setting **batch** in the metadata to the same batch.
2. Set **phenotype** of the CRISPR-edited sample to **tumor**.
3. Set the **phenotype** of the non-edited sample to **normal**.

And kick off the **variant2** pipeline, the normal whole genome sequencing pipeline. An example YAML template is below:

```yaml
---
details:
  - analysis: variant2
    genome_build: hg38
    algorithm:
      aligner: bwa
      variantcaller: mutect2
      tools_on: [gemini]
```

And an example metadata file:

```csv
samplename,description,batch,phenotype,sex,cas9,gRNA
Hs27_HSV1.cram,Hs27_HSV1,noCas9_nogRNA,normal,male,no,yes
Hs27_HSV1_Cas9.cram,Hs27_HSV1_Cas9,noCas9,normal,male,yes,no
Hs27_HSV1_UL30_5.cram,Hs27_HSV1_UL30_5,noCas9_nogRNA,tumor,male,yes,yes
Hs27_HSV1_UL30_5_repeat.cram,Hs27_HSV1_UL30_5_repeat,noCas9,tumor,male,yes,yes
```

# Find predicted off-target sites
There are several tools to do this, a common one folks use is cas-offinder, so
that is what we will use. There is a [web app](http://www.rgenome.net/cas-offinder/) but it will only return 1,000 events per
class. Usually this is fine, but if you allow bulges you can get a lot more offtarget sites so you might bump into this limit.

First install cas-offinder, there is a conda package so this is easy:

```bash
conda create -n crispr -c bioconda cas-offinder
```

There is a companion python wrapper cas-offinder-bulge that can also predict
offtarget sites taking bulges into effect. You can download it
[here](https://raw.githubusercontent.com/hyugel/cas-offinder-bulge/master/cas-offinder-bulge) if
you need to do that.

You'll need to know the sequence of one or more guides you want to check. You will also need to know
the PAM sequence for the endonuclease that is being used. 

You can run cas-offinder like this:

```bash
cas-offinder input.txt C output.txt
```

where input.txt has this format:

```
hg38.fa
NNNNNNNNNNNNNNNNNNNNNNNGRRT
ACACGTGAAAGACGGTGACGGNNGRRT 6
```

`hg38.fa` is the path to the FASTA file of the hg38 genome. NNNNNNNNNNNNNNNNNNNNNNNNGRRT is the length of the guide sequence you are interested in with the PAM sequence tacked on the end.
ACACACGTGAAAGACGGTGACGGNNGRRT is the guide sequence with the PAM sequence tacked on the end. 6 is the number of mismatches you are allowing here, it will look for sites with that many
or less mismatches.

If you want to look for bulges, use cas-offinder-bulge with this format:

```
hg38.fa
NNNNNNNNNNNNNNNNNNNNNNNGRRT 2 1
ACACGTGAAAGACGGTGACGGNNGRRT 6
```

Where the 2 says to look for a DNA bulge and the 1 a RNA bulge. You can do one or the other, neither or both.

After you run cas-offinder you can conver the output to a sorted BED file for use with intersecting the your variants:

```bash
cat output.txt | sed 1d | awk '{printf("%s\t%s\t%s\n",$4, $5-10,$5+10)}' | sort -V -k 1,1 -k2,2n  > output.bed
```

# Overlap variants
Finally, use the BED file of predicted off-target sites to pull out possible off-target variant calls:

```bash
bedtools intersect -header -u -a noCas9_nogRNA-mutect2-annotated.vcf.gz -b output.bed
```

And you are done!
