8_7_25 Billingsley

### Quick Start guide to aligning a bulk RNA-Seq experiment using Nextflow pipeline on fas.rc

For review (from Alex Bartlett)

https://hbc.github.io/Platform/<br>
https://hbc.github.io/Platform/pipelines/#nextflow-in-fas<br>
https://nf-co.re/rnaseq/3.19.0<br><br><br><br>

Make a directory on fas.rc on netscratch, like this example:

/n/netscratch/hsph_bioinfo/Lab/sarah_hill/<br><br><br><br>



Make subdirectories for the fastqs, and the nextflow run, like this:

/n/netscratch/hsph_bioinfo/Lab/sarah_hill/250701_SH13270_fastq/

/n/netscratch/hsph_bioinfo/Lab/sarah_hill/nfcore_run_jb <br><br><br><br>



From fas.rc, copy your fastqs from o2 to the fastq directory like this:

rsync -chavzP jmb17@Transfer.RC.hms.harvard.edu:/n/data1/cores/bcbio/PIs/sarah_hill/hill_rnaseq_human_endometrial_cancer_hbc05261/data/ /n/netscratch/hsph_bioinfo/Lab/sarah_hill/250701_SH13270_fastq/<br><br><br><br>


For a Nextflow run you'll need 5 subdirectories in your nfcore_run directory. (examples from ZhuZhou and UpenBhattarai)

drwxr-sr-x. 2 zzhuo hsph_bioinfo 4096 Jun 27 17:52 config<br>
drwxr-sr-x. 2 zzhuo hsph_bioinfo 4096 Jun 27 17:42 data<br>
drwxr-sr-x. 2 zzhuo hsph_bioinfo 4096 Jul  3 17:53 nfcore_output<br>
drwxr-sr-x. 2 zzhuo hsph_bioinfo 4096 Jun 27 17:38 resources<br>
drwxr-sr-x. 2 zzhuo hsph_bioinfo 4096 Jul  3 17:53 workdir<br>
<br><br>

in /config, you'll need 4 files

-rw-r-xr--. 1 zzhuo hsph_bioinfo   90 Jun 27 17:37 fas.config<br>
-rw-r--r--. 1 zzhuo hsph_bioinfo 3276 Jun 27 17:52 rnaseq.json<br>
-rw-r--r--. 1 zzhuo hsph_bioinfo 4007 Jun 27 17:37 rnaseq.resources.config<br>
-rw-r--r--. 1 zzhuo hsph_bioinfo  840 Jun 27 17:50 rrna-db-defaults.txt<br>
<br><br>

fas.config is this:

executor.$slurm.queueSize = 30
process {
    executor = 'slurm'
    queue    = 'shared'
}

<br><br>

rrna-db-defaults.txt is this:


https://raw.githubusercontent.com/biocore/sortmerna/v4.2.0/data/rRNA_databases/rfam-5.8s-database-id98.fasta
https://raw.githubusercontent.com/biocore/sortmerna/v4.2.0/data/rRNA_databases/rfam-5s-database-id98.fasta
https://raw.githubusercontent.com/biocore/sortmerna/v4.2.0/data/rRNA_databases/silva-arc-16s-id95.fasta
https://raw.githubusercontent.com/biocore/sortmerna/v4.2.0/data/rRNA_databases/silva-arc-23s-id98.fasta
https://raw.githubusercontent.com/biocore/sortmerna/v4.2.0/data/rRNA_databases/silva-bac-16s-id90.fasta
https://raw.githubusercontent.com/biocore/sortmerna/v4.2.0/data/rRNA_databases/silva-bac-23s-id98.fasta
https://raw.githubusercontent.com/biocore/sortmerna/v4.2.0/data/rRNA_databases/silva-euk-18s-id95.fasta
https://raw.githubusercontent.com/biocore/sortmerna/v4.2.0/data/rRNA_databases/silva-euk-28s-id98.fasta

<br><br>
rnaseq.json: see below

<br><br>
rnaseq.resources.config: see below

<br><br><br><br>


in /data, place a metadata.csv file. This points to the fastqs, and you can include whatever metadata you want for each sample.

A minimal example requires 4 columns, sample, fastq_1, fastq_2, strandedness

A test run on a couple of samples looks like this for my run:


sample,fastq_1,fastq_2,strandedness
Patient_10,/n/netscratch/hsph_bioinfo/Lab/sarah_hill/250701_SH13270_fastq/20250701_Patient_10_SH13270_S351_L008_R1_001.fastq.gz,/n/netscratch/hsph_bioinfo/Lab/sarah_hill/250701_SH13270_fastq/20250701_Patient_10_SH13270_S351_L008_R2_001.fastq.gz,auto
Patient_11,/n/netscratch/hsph_bioinfo/Lab/sarah_hill/250701_SH13270_fastq/20250701_Patient_11_SH13270_S363_L008_R1_001.fastq.gz,/n/netscratch/hsph_bioinfo/Lab/sarah_hill/250701_SH13270_fastq/20250701_Patient_11_SH13270_S363_L008_R2_001.fastq.gz,auto
<br><br>

Here are the columns for a run Zhu did, lots of metadata!: sample,fastq_1,fastq_2,strandedness,phase,Donor_id,Subject_code,Sex,Age,Number_of_pills,Biopsies_collection_date,Biopsies_collection_time,RNA_extraction_date,RNA_conc(ng/ul),A260/A280
<br><br><br><br>



In /resources, place the .fa and .gtf files, for my run,  gencode.v47.annotation.gtf and GRCh38.p14.genome.fa
<br><br><br><br>


In /workdir, place your slurm script.  

In /workdir I also placed a /work subdirectory, I suspect this would be created and populated by the run if you don't create it yourself, I'm not sure.
<br><br><br><br>



My slurm script was this, based on scripts from Zhu and Upen, I just changed some paths:

```#!/bin/bash

#SBATCH --job-name=Nextflow      # Job name
#SBATCH --partition=shared            # Partition name
#SBATCH --time=3-00:00                 # Runtime in D-HH:MM format
#SBATCH --nodes=1                      # Number of nodes (keep at 1)
#SBATCH --ntasks=1                     # Number of tasks per node (keep at 1)
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G                     # Memory needed per node (total)
#SBATCH --error=%x_%j.err           # File to which STDERR will be written, including job ID
#SBATCH --output=%x_%j.out          # File to which STDOUT will be written, including job ID
#SBATCH --mail-type=ALL                # Type of email notification (BEGIN, END, FAIL, ALL)

module load jdk/21.0.2-fasrc01

export NXF_APPTAINER_CACHEDIR=/n/netscratch/hsph_bioinfo/Lab/sarah_hill/nfcore_run_jb/workdir/tmp
export NXF_SINGULARITY_LIBRARYDIR=/n/holylfs05/LABS/hsph_bioinfo/Lab/shared_resources/nextflow/nfcore-rnaseq


/n/holylfs05/LABS/hsph_bioinfo/Lab/shared_resources/nextflow/nextflow run nf-core/rnaseq -r 3.14.0 \
  -profile singularity \
  -params-file /n/netscratch/hsph_bioinfo/Lab/sarah_hill/nfcore_run_jb/config/rnaseq.json \
  -c /n/netscratch/hsph_bioinfo/Lab/sarah_hill/nfcore_run_jb/config/fas.config \
  -c /n/netscratch/hsph_bioinfo/Lab/sarah_hill/nfcore_run_jb/config/rnaseq.resources.config \
  --input /n/netscratch/hsph_bioinfo/Lab/sarah_hill/nfcore_run_jb/data/metadata.csv \
  --outdir /n/netscratch/hsph_bioinfo/Lab/sarah_hill/nfcore_run_jb/nfcore_output \
  -resume
  ```
<br><br><br><br>




Here's the rnsseq.json file to place in /config. I copied this from Zhu, and just changed three file paths in this to point at the files in my own /resources and /config directories,  nothing else was changed.



```
{
    "remove_ribo_rna": false,
    "custom_config_base": "https://raw.githubusercontent.com/nf-core/configs/master",
    "skip_deseq2_qc": false,
    "gencode": true,
    "monochrome-logs": false,
    "validation-S3Path-check": false,
    "umitools_dedup_stats": false,
    "plaintext_email": false,
    "save_reference": true,
    "skip_markduplicates": true,
    "ribo_database_manifest": "/n/netscratch/hsph_bioinfo/Lab/sarah_hill/nfcore_run_jb/config/rrna-db-defaults.txt",
    "monochrome_logs": false,
    "aligner": "star_salmon",
    "max_cpus": 16,
    "featurecounts_group_type": "gene_biotype",
    "save_bbsplit_reads": false,
    "skip_multiqc": false,
    "validationFailUnrecognisedParams": false,
    "skip_preseq": true,
    "skip_dupradar": true,
    "save_align_intermeds": false,
    "gtf": "/n/netscratch/hsph_bioinfo/Lab/sarah_hill/nfcore_run_jb/resources/gencode.v47.annotation.gtf",
    "max_multiqc_email_size": "25.MB",
    "max_time": "240.h",
    "save_trimmed": false,
    "salmon_quant_libtype": "A",
    "pseudo_aligner": "salmon",
    "min_trimmed_reads": 10000,
    "skip_fastqc": false,
    "pseudo_aligner_kmer_size": 31,
    "deseq2_vst": true,
    "umitools_extract_method": "string",
    "validate_params": true,
    "validationShowHiddenParams": false,
    "bam_csi_index": false,
    "with_umi": false,
    "skip_qc": false,
    "version": false,
    "publish_dir_mode": "copy",
    "validationS3PathCheck": false,
    "input": "",
    "skip_bigwig": true,
    "igenomes_base": "s3://ngi-igenomes/igenomes",
    "validationSkipDuplicateCheck": false,
    "validation-show-hidden-params": false,
    "validationSchemaIgnoreParams": "genomes,igenomes_base",
    "validation-lenient-mode": false,
    "skip_gtf_transcript_filter": false,
    "stringtie_ignore_gtf": false,
    "skip_umi_extract": true,
    "kallisto_quant_fraglen_sd": 200,
    "validation-skip-duplicate-check": false,
    "save_unaligned": false,
    "featurecounts_feature_type": "exon",
    "skip_gtf_filter": false,
    "custom_config_version": "master",
    "fasta": "/n/netscratch/hsph_bioinfo/Lab/sarah_hill/nfcore_run_jb/resources/GRCh38.p14.genome.fa",
    "max_memory": "180.GB",
    "hisat2_build_memory": "250.GB",
    "skip_rseqc": true,
    "save_non_ribo_reads": false,
    "monochromeLogs": false,
    "skip_alignment": false,
    "skip_qualimap": false,
    "star_ignore_sjdbgtf": false,
    "skip_stringtie": true,
    "gtf_extra_attributes": "gene_name",
    "save_merged_fastq": false,
    "save_umi_intermeds": false,
    "min_mapped_reads": 5,
    "rseqc_modules": "bam_stat,inner_distance,infer_experiment,junction_annotation,junction_saturation,read_distribution,read_duplication",
    "validationLenientMode": false,
    "gtf_group_features": "gene_id",
    "skip_trimming": false,
    "kallisto_quant_fraglen": 200,
    "igenomes_ignore": false,
    "skip_pseudo_alignment": true,
    "outdir": "REPLACE_ME",
    "skip_bbsplit": true,
    "help": false,
    "validation-fail-unrecognised-params": false,
    "test_data_base": "https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3",
    "validation-schema-ignore-params": "genomes,igenomes_base",
    "umitools_grouping_method": "directional",
    "skip_biotype_qc": false
}
    
```
<br><br><br><br>

Here's the rnaseq.resources.config file, copied from Zhu, nothing was changed.


  
  ```
process {
  withName: 'NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:BAM_SORT_STATS_SAMTOOLS:BAM_STATS_SAMTOOLS:SAMTOOLS_FLAGSTAT' {
    cpus = { 3 * task.attempt }
    memory = { 5.GB * task.attempt }
  }
  withName: 'NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:BAM_SORT_STATS_SAMTOOLS:BAM_STATS_SAMTOOLS:SAMTOOLS_IDXSTATS' {
    cpus = { 1 * task.attempt }
    memory = { 5.GB * task.attempt }
  }
  withName: 'NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:BAM_SORT_STATS_SAMTOOLS:BAM_STATS_SAMTOOLS:SAMTOOLS_STATS' {
    cpus = { 2 * task.attempt }
    memory = { 5.GB * task.attempt }
  }
  withName: 'NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:BAM_SORT_STATS_SAMTOOLS:SAMTOOLS_INDEX' {
    cpus = { 2 * task.attempt }
    memory = { 5.GB * task.attempt }
  }
  withName: 'NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:BAM_SORT_STATS_SAMTOOLS:SAMTOOLS_SORT' {
    cpus = { 5 * task.attempt }
    memory = { 10.GB * task.attempt }
  }
  withName: 'NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN' {
    cpus = { 5 * task.attempt }
    memory = { 40.GB * task.attempt }
  }
  withName: 'NFCORE_RNASEQ:RNASEQ:CUSTOM_DUMPSOFTWAREVERSIONS' {
    cpus = { 1 * task.attempt }
    memory = { 5.GB * task.attempt }
  }
  withName: 'NFCORE_RNASEQ:RNASEQ:DESEQ2_QC_STAR_SALMON' {
    cpus = { 2 * task.attempt }
    memory = { 5.GB * task.attempt }
  }
  withName: 'NFCORE_RNASEQ:RNASEQ:FASTQ_FASTQC_UMITOOLS_TRIMGALORE:FASTQC' {
    cpus = { 2 * task.attempt }
    memory = { 5.GB * task.attempt }
  }
  withName: 'NFCORE_RNASEQ:RNASEQ:FASTQ_FASTQC_UMITOOLS_TRIMGALORE:TRIMGALORE' {
    cpus = { 4 * task.attempt }
    memory = { 20.GB * task.attempt }
  }
  withName: 'NFCORE_RNASEQ:RNASEQ:FASTQ_SUBSAMPLE_FQ_SALMON:FQ_SUBSAMPLE' {
    cpus = { 1 * task.attempt }
    memory = { 5.GB * task.attempt }
  }
  withName: 'NFCORE_RNASEQ:RNASEQ:FASTQ_SUBSAMPLE_FQ_SALMON:SALMON_QUANT' {
    cpus = { 2 * task.attempt }
    memory = { 40.GB * task.attempt }
  }
  withName: 'NFCORE_RNASEQ:RNASEQ:MULTIQC' {
    cpus = { 2 * task.attempt }
    memory = { 5.GB * task.attempt }
  }
  withName: 'NFCORE_RNASEQ:RNASEQ:MULTIQC_CUSTOM_BIOTYPE' {
    cpus = { 1 * task.attempt }
    memory = { 5.GB * task.attempt }
  }
  withName: 'NFCORE_RNASEQ:RNASEQ:PREPARE_GENOME:CUSTOM_GETCHROMSIZES' {
    cpus = { 1 * task.attempt }
    memory = { 5.GB * task.attempt }
  }
  withName: 'NFCORE_RNASEQ:RNASEQ:PREPARE_GENOME:GTF2BED' {
    cpus = { 1 * task.attempt }
    memory = { 16.GB * task.attempt }
  }
  withName: 'NFCORE_RNASEQ:RNASEQ:PREPARE_GENOME:GTF_FILTER' {
    cpus = { 1 * task.attempt }
    memory = { 5.GB * task.attempt }
  }
  withName: 'NFCORE_RNASEQ:RNASEQ:PREPARE_GENOME:MAKE_TRANSCRIPTS_FASTA' {
    cpus = { 1 * task.attempt }
    memory = { 5.GB * task.attempt }
  }
  withName: 'NFCORE_RNASEQ:RNASEQ:QUALIMAP_RNASEQ' {
    cpus = { 1 * task.attempt }
    memory = { 20.GB * task.attempt }
  }
  withName: 'NFCORE_RNASEQ:RNASEQ:QUANTIFY_STAR_SALMON:SALMON_QUANT' {
    cpus = { 5 * task.attempt }
    memory = { 30.GB * task.attempt }
  }
  withName: 'NFCORE_RNASEQ:RNASEQ:QUANTIFY_STAR_SALMON:SE_GENE' {
    cpus = { 2 * task.attempt }
    memory = { 5.GB * task.attempt }
  }
  withName: 'NFCORE_RNASEQ:RNASEQ:QUANTIFY_STAR_SALMON:SE_GENE_LENGTH_SCALED' {
    cpus = { 2 * task.attempt }
    memory = { 5.GB * task.attempt }
  }
  withName: 'NFCORE_RNASEQ:RNASEQ:QUANTIFY_STAR_SALMON:SE_GENE_SCALED' {
    cpus = { 3 * task.attempt }
    memory = { 5.GB * task.attempt }
  }
  withName: 'NFCORE_RNASEQ:RNASEQ:QUANTIFY_STAR_SALMON:SE_TRANSCRIPT' {
    cpus = { 2 * task.attempt }
    memory = { 5.GB * task.attempt }
  }
  withName: 'NFCORE_RNASEQ:RNASEQ:QUANTIFY_STAR_SALMON:TX2GENE' {
    cpus = { 1 * task.attempt }
    memory = { 5.GB * task.attempt }
  }
  withName: 'NFCORE_RNASEQ:RNASEQ:QUANTIFY_STAR_SALMON:TXIMPORT' {
    cpus = { 2 * task.attempt }
    memory = { 6.GB * task.attempt }
  }
  withName: 'NFCORE_RNASEQ:RNASEQ:SUBREAD_FEATURECOUNTS' {
    cpus = { 3 * task.attempt }
    memory = { 5.GB * task.attempt }
  }
  errorStrategy = 'retry'
  maxRetries = 2
}
  
 ```





