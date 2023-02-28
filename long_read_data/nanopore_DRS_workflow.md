# Workflow for processing Nanopore Direct RNA sequencing

_Zhu Zhuo_ \
_Feb 28, 2023_

This document briefly describes the steps and tips for processing the Nanopore Direct RNA Sequencing data on HMS O2 cluster.

It includes basecalling, filtering low quality reads, mapping to genome and transcriptome references, checking alignment stats, and getting read counts.

Feel free to contact me if there is any question.

## 1. **Tidy up fast5 files**

 Create a folder called fast5 and then create a folder for each sample that contains all the fast5 files within fast5.

## 2. Run Guppy-gpu for basecalling.

Guppy can't be installed without admin. So I used a Guppy-gpu containers for basecalling.

### 2a. Import Guppy containers in O2

O2 has Singularity installed and Singularity is fully compatible with existing Docker images. Singularity can only be used to run images that have been tested and approved. The testing process is fully automated, and can be initiated by any users. I used the Guppy-gpu docker repository with most downloads (i.e. genomicpariscentre/guppy-gpu), and followed the intruction here to import Guppy-gpu: https://harvardmed.atlassian.net/wiki/spaces/O2/pages/1616511103/Running+Singularity+Containers+in+O2

All system libraries inside the Docker container needs to be the last version available, otherwise the Docker container will fail to import. Usually this can be done by running the command apt update and apt-get upgrade at the end of the installation process. I tried but it doesn't work. I had to do the following to get it work (can't do this without Victor's help):

```
## locally in a termial
docker pull genomicpariscentre/guppy-gpu
docker image ls
docker run -it genomicpariscentre/guppy-gpu
## then in the container
apt-get update -y
apt-get upgrade -y

## locally, open another terminal, push the updated container to my personal docker account
docker container ls
docker commit 12051984a63c zhuzhuo/guppy-gpu_updated
docker image push zhuzhuo/guppy-gpu_updated

## on o2, submit the container
csubmitter --name guppy-gpu_updated --image-uri docker://zhuzhuo/guppy-gpu_updated:latest ## Passed!
```

### 2b. Running Guppy-gpu for basecalling

It's worthwhile to specify in `--gres` to use teslaV100 for the job: `#SBATCH --gres=gpu:teslaV100:2`

teslaK80 might not be working https://gist.github.com/sirselim/13f70ae69f2a512e7d9e1f00f9704f53. It's possible that's why the job array failed when I did `#SBATCH --gres=gpu:2`, because some type of GPU is not working.

Use `nvidia-smi` to check what GPU is used.

Below is an example of the script. I disabled q score filtering.

```
#!/bin/bash
#SBATCH -c 4
#SBATCH -t 12:00:00
#SBATCH -p gpu
#SBATCH --gres=gpu:teslaV100:2
#SBATCH -J guppy
#SBATCH -o %x_%j.o
#SBATCH -e %x_%j.e

module load gcc/9.2.0 cuda/11.2

singularity exec --nv /n/app/singularity/containers/hmsid/guppy-gpu_updated.sif guppy_basecaller -i /path/to/fast5/sample -s //path/to/fastq/sample -c rna_r9.4.1_70bps_hac.cfg --compress_fastq --disable_qscore_filtering --records_per_fastq 0 -x "cuda:0,1"
```

## 3. step 3 - all steps after guppy basecalling

Below is a job array script for all the steps after basecalling, but it can be modified to accommodate your needs.

Multiple software are used in this step. I installed some of the software in a conda environment I created (called conda-env-python3) and some of them in my personal software folder.

Note: .ipynb_checkpoints folder need to be removed or `source activate` is not working!

Step 3 includes the following substeps:

a. Using Nanoplot and MinIONQC to QC the fastq data.

b. Filtering out reads with low quality scores using NanoFilt. Cutoff of Q score = 5 is used in the script below.

c. Mapping the reads to the genome and the transcriptome using Minimap2.

d. Using NanoPlot to check the bam files.

I also use this commands below to check mapped percentage and mapping identity:

```
for f in *.bam; do if [ -s $f ]; then echo $f `samtools stats -F3840 $f | awk '{if ($2=="sequences:"){nseq=$3} else if($2=="reads" && $3=="mapped:"){printf("mapped: %s (%.1f%)\n",$4,100*$4/nseq)} else if($2=="error" && $3=="rate:"){printf("identity: %.2f%\n",100-$4*1e2)}}'`; fi; done
```

e. Get the read counts from the alignments to genome and transcriptome using featureCounts and Nanocount, respectively.

We also did a test with `--primary` turned on for featureCounts. The default of featureCounts only counts uniquely mapped reads by checking NH:i tag, but minimap2 doesn't output this tag. With `--primary` on, featureCounts count each read only once.


```
#!/bin/bash
#SBATCH -J step3
#SBATCH -p short
#SBATCH -c 12
#SBATCH --mem=128G
#SBATCH -t 2:00:00
#SBATCH -o %x_%j.o
#SBATCH -e %x_%j.e
#SBATCH --array=1-3%2

base=/path/to/base

samples=(s1 s2 s3)
sample=${samples[$SLURM_ARRAY_TASK_ID-1]}

echo "#QC of the fastq#"

echo "##Nanoplot of the fastq starting##"

source /home/${USER}/.bashrc
source activate conda-env-python3

NanoPlot -t 12 --summary ${base}/fastq/${sample}/sequencing_summary.txt --loglength -o ${base}/fastq/nanoplot/${sample}

echo "##Nanoplot of the fastq done##"

echo "##MinIONQC of the fastq starting##"

module load gcc/9.2.0 R/4.2.1

Rscript ~/bin/MinIONQC.R -p 12 -q 5 -i ${base}/fastq/${sample}/sequencing_summary.txt -o ${base}/fastq/minIONQC/${sample}

echo "##MinIONQC of the fastq done##"


echo "#Filter fastqs#"

echo "##NanoFilt starting##"

cat ${base}/fastq/${sample}/fastq_*.fastq.gz | gunzip -c | NanoFilt -q 5 -s ${base}/fastq/${sample}/sequencing_summary.txt --readtype 1D | gzip > ${base}/fastq_filt/${sample}.filt.fq.gz

echo "##NanoFilt done##"


echo "#Mapping#"

echo "##Mapping to genome reference starting##"

export PATH=~/bin/minimap2:$PATH

minimap2 -ax splice -uf -k14 -t 12 ${base}/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa ${base}/fastq_filt/${sample}.filt.fq.gz | samtools sort -O BAM -@12 -o ${base}/aln/${sample}.g.bam -

samtools index ${base}/aln/${sample}.g.bam

echo "##Mapping to genome reference done##"

echo "##Mapping to transcriptome reference starting##"

export PATH=~/bin/minimap2:$PATH

minimap2 -ax map-ont -N 10 -t 12 ${base}/ref/transcripts_union.fa.gz ${base}/fastq_filt/${sample}.filt.fq.gz | samtools sort -O BAM -@12 -o ${base}/aln/${sample}.t.bam -

samtools index ${base}/aln/${sample}.t.bam

echo "##Mapping to transcriptome reference done##"

echo "##Nanoplot of the bam starting##"

NanoPlot -t 12 --bam ${base}/aln/${sample}.g.bam -o ${base}/aln/nanoplot/${sample}.g
NanoPlot -t 12 --bam ${base}/aln/${sample}.t.bam -o ${base}/aln/nanoplot/${sample}.t

echo "##Nanoplot of the bam done##"


echo "#Getting read counts#"

featureCounts -T 12 -L -a ${base}/ref/Homo_sapiens.GRCh38.107.gtf.gz -o ${base}/counts/gx_counts/${sample}.featureCount.txt ${base}/aln/${sample}.g.bam
NanoCount -i ${base}/aln/${sample}.t.bam -o ${base}/counts/tx_counts/${sample}.tx.tsv

echo "#Getting read counts done#"
```
