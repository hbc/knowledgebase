---
tags:
title: Genome Assembly Using PacBio Reads Only
author: Zhu Zhuo
created: '2019-09-13'
---

# Genome Assembly Using PacBio Reads Only

This tutorial is based on a bacterial genome assembly project using PacBio sequencing reads only, but it can be followed for genome assembly of other species or using other type long reads with no or a little modification.

## Demultiplex

If the sequencing core hasn't demultiplexed the data, [`lima`](https://github.com/PacificBiosciences/barcoding) can be used for demultiplexing.

## Convert `bam` file to `fastq` file

`subreads.bam` files contain the subreads data and we will convert it from `bam` format to `fastq` format as most assemblers take `fastq` as input.

[A note on the output from PacBio:](https://pacbiofileformats.readthedocs.io/en/5.1/Primer.html)
> Unaligned BAM files representing the subreads will be produced natively by the PacBio instrument. The subreads BAM will be the starting point for secondary analysis. In addition, the scraps arising from cutting out adapter and barcode sequences will be retained in a `scraps.bam` file, to enable reconstruction of HQ regions of the ZMW reads, in case the customer needs to rerun barcode finding with a different option.

Below is an example of slurm script using `bedtools bamtofastq` to convert `bam` to `fastq`.

```
#!/bin/sh
#SBATCH -p medium
#SBATCH -J bam2fq
#SBATCH -o x%_%j.o
#SBATCH -e x%_%j.e
#SBATCH -t 00-23:59:00
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH --array=1-n%5 #change n to the number of subreads.bam files

module load bedtools/2.27.1

files=(/path/to/subreads.bam)
file=${files[$SLURM_ARRAY_TASK_ID-1]}
sample=`basename file .bam`

echo $file
echo $sample

bedtools bamtofastq -i $file -fq $sample".fq"
```

PacBio Sequel Sequencer reports all base qualities as PHRED 0 (ASCII !). So the quality score for sequel data are all `!` in the `fastq` file generated.

## Genome assembly

### Using Canu for genome assembly

[Canu](https://github.com/marbl/canu) does correction, trimming and assembly in a single command. Follow its github page to install the software to a location of preference.

An example of slurm script for a single sample:
```
#!/bin/sh
#SBATCH -p priority
#SBATCH -J canu
#SBATCH -o %x_%j.o
#SBATCH -e %x_%j.e
#SBATCH -t 0-23:59:00
#SBATCH -c 1
#SBATCH --mem=1G

module load java/jdk-1.8u112

export PATH=/path/to/canu:$PATH

canu -p sampleName -d sampleName genomeSize=7m \
  gridOptions="--time=1-23:59:00 --partition=medium" \
  -pacbio-raw /path/to/converted.fq
```
This is for the 'master' job, so only 1 CPU and 1 Gb memory should be sufficient. Canu will evaluate the resources available and automatically submit jobs to the queue in `gridOptions`.

_Note_: An alternative is to install a bioconda recipe. But the conda verison is not up-to-date and some additional parameters may need to be specified in the command.

### Using Unicycler for bacterial genome assembly

Follow [Unicycler](https://github.com/rrwick/Unicycler#method-long-read-only-assembly) instructions. [Racon](https://github.com/isovic/racon) is also required and should be installed before running Unicycler.

```
module load gcc/6.2.0 python/3.6.0
git clone https://github.com/rrwick/Unicycler.git
cd Unicycler
python3 setup.py install --user
```
Example of a slurm script
```
#!/bin/sh
#SBATCH -p priority
#SBATCH -J unicycler
#SBATCH -o %x_%j.o
#SBATCH -e %x_%j.e
#SBATCH -t 29-23:59:00
#SBATCH -c 20
#SBATCH --mem=200G

module load gcc/6.2.0 python/3.6.0 bowtie2/2.3.4.3 samtools/1.9 blast/2.6.0+
export PATH=/path/to/racon:$PATH

/path/to/unicycler-runner.py -l /path/to/fastq -t 20 -o sample
```
## Assembly quality

### Basic assembly metrics

Download [Quast](http://bioinf.spbau.ru/quast) for basic assembly metrics, such as total length, number of contigs and N50.
`/path/to/quast.py -o output_folder -t 6 assembly.fa`

### Assembly completeness

Use [BUSCO](https://busco.ezlab.org/) for evaluating the completeness of the genome assembly. BUSCO has a lot of dependencies, so it is better to install a conda recipe.

```
source activate conda-env # activate conda environment. If you don't have one, you may need to create a conda environment.
conda install -c bioconda busco
conda deactivate # deactivate conda environment
```
Download the BUSCO database for the species and run BUSCO

`run_busco -i assembly.fa -o output_folder -l species_odb -m geno`

BUSCO is also the abbreviation for Benchmarking Universal Single-Copy Orthologs, which is single-copy ortholog found in >90% of species. The more BUSCOs are present, the more complete the genome assembly is.
