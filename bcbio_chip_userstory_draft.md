# ChIP-seq

bcbio evaluates data quality using FASTQC and filters and trims reads as necessary. Reads are aligned to the genome and only uniquely aligning reads are retained. Peaks are called with MACS2. A consensus peak file is generated using bedops and this consensus peak list is used to generate a count matrix using featureCounts.

## Description of example dataset

We will be working using ChIP-seq data from a recent publication in Neuron by Baizabal et al. (2018) [[1]](https://doi.org/10.1016/j.neuron.2018.04.033). 
Baizabal et al. sought to understand how chromatin-modifying enzymes function in neural stem cells to establish the epigenetic landscape that determines cell type and stage-specific gene expression. Chromatin-modifying enzymes are transcriptional regulators that control gene expression through covalent modification of DNA or histones. 

This specific dataset focuses on the transcriptional regulator **PRDM16, which is a chromatin-modifying enzyme** that belongs to the larger PRDM (Positive Regulatory Domain) protein family, that is structurally defined by the **presence of a conserved N-terminal histone methyltransferase PR domain** ([Hohenauer and Moore, 2012](https://journals.biologists.com/dev/article/139/13/2267/45169/The-Prdm-family-expanding-roles-in-stem-cells-and)). The authors generated ChIP-seq data for PRDM16. Our dataset consists of two WT samples and two KO samples.


### 1. Download the example data and configuration files
This downloads the input data, creates the project structure and example configuration files.

#### 1.1 Create input directory and download FASTQ files.

To download the files used in this experiment, you will need to use the SRA Toolkit. Once installed, you can use the following commands:

```bash
mkdir chip-example
cd chip-example
mkdir -p fastq

# For WT Sample 1 ChIP
fastq-dump --accession SRR6823762
fastq-dump --accession SRR6823763
cat SRR6823762.fastq SRR6823763.fastq | gzip > wt_sample1_chip.fastq.gz
rm SRR6823762.fastq SRR6823763.fastq

# For WT Sample 1 Input
fastq-dump --accession SRR6823764
fastq-dump --accession SRR6823765
cat SRR6823764.fastq SRR6823765.fastq | gzip > wt_sample1_input.fastq.gz
rm SRR6823764.fastq SRR6823765.fastq

# For WT Sample 2 ChIP
fastq-dump --accession SRR6823766
fastq-dump --accession SRR6823767
cat SRR6823766.fastq SRR6823767.fastq | gzip > wt_sample2_chip.fastq.gz
rm SRR6823766.fastq SRR6823767.fastq

# For WT Sample 2 Input
fastq-dump --accession SRR6823768
fastq-dump --accession SRR6823769
cat SRR6823768.fastq SRR6823769.fastq | gzip > wt_sample2_input.fastq.gz
rm SRR6823768.fastq SRR6823769.fastq

# For KO Sample 1 ChIP
fastq-dump --accession SRR6823770
cat SRR6823770.fastq | gzip > ko_sample1_chip.fastq.gz
rm SRR6823770.fastq

# For KO Sample 1 Input
fastq-dump --accession SRR6823771
cat SRR6823771.fastq | gzip > ko_sample1_input.fastq.gz
rm SRR6823771.fastq

# For KO Sample 2 ChIP                                                                    
fastq-dump --accession SRR6823772
cat SRR6823772.fastq | gzip > ko_sample2_chip.fastq.gz
rm SRR6823772.fastq

# For KO Sample 2 Input
fastq-dump --accession SRR6823773
cat SRR6823773.fastq | gzip > ko_sample2_input.fastq.gz
rm SRR6823773.fastq
```

#### 1.2 Create the template YAML file describing the ChIP-seq analysis

Create the yaml file and then copy in the example text.

```bash
mkdir templates
cd templates
nano chip-example.yaml
```

chip-example.yaml: 

```yaml
details:
- analysis: chip-seq
  genome_build: mm10
  algorithm:
    aligner: bowtie2
    adapters: illumina
    trim_reads: read_through
    peakcaller: macs2
    chip_method: chip
    keep_duplicates: False
    keep_multimapped: False
upload:
  dir: ../final
```

Here we are telling bcbio that we want our reads trimmed (trim_reads: read_through) rather then letting bcbio decide.

#### 1.3 Create a sample sheet

Create the sample sheet then copy in the example text.

```bash
mkdir metadata
cd metadata
nano neurons.csv
```

```
File,description,genotype,enzyme,batch,phenotype,antibody
ko_sample1_chip.fastq.gz,ko_sample1_chip,KO,PRDM19,pair1,chip,narrow
ko_sample1_input.fastq.gz,ko_sample1_input,KO,INPUT,pair1,input,narrow
ko_sample2_chip.fastq.gz,ko_sample2_chip,KO,PRDM19,pair2,chip,narrow
ko_sample2_input.fastq.gz,ko_sample2_input,KO,INPUT,pair2,input,narrow
wt_sample1_chip.fastq.gz,wt_sample1_chip,WT,PRDM19,pair3,chip,narrow
wt_sample1_input.fastq.gz,wt_sample1_input,WT,INPUT,pair3,input,narrow
wt_sample2_chip.fastq.gz,wt_sample2_chip,WT,PRDM19,pair4,chip,narrow
wt_sample2_input.fastq.gz,wt_sample2_input,WT,INPUT,pair4,input,narrow
```

The necessary columns here are: `File`, `description`, `batch`, `phenotype` and `antibody`. 
 `batch`,`phenotype`, and `antibody` are columns unique to ChIP-seq.

`batch` matches your input samples with their respective chips and the `phenotype` column tells bcbio if a sample is an input or chip.

Here we have one input for every chip. For example ko_sample1 has ko_sample1_chip and ko_sample1_input. These are pair1. However, sometimes the same input is used for multiple chips. Here is the same file but assuming that we also ran a h3k4me1 chip on all samples:


```
File,description,genotype,enzyme,batch,phenotype,antibody
ko_sample1_chip_prdm.fastq.gz,ko_sample1_chip_prdm,KO,PRDM19,pair1,chip,narrow
ko_sample1_chip_h3k4.fastq.gz,ko_sample1_chip_h3k4,KO,H3K4Me1,pair2,chip,narrow
ko_sample1_input.fastq.gz,ko_sample1_input,KO,INPUT,pair1;pair2,input,narrow
ko_sample2_chip_prdm.fastq.gz,ko_sample2_chip_prdm,KO,PRDM19,pair3,chip,narrow
ko_sample2_chip_h3k4.fastq.gz,ko_sample2_chip_h3k4,KO,H3K4Me1,pair4,chip,narrow
ko_sample2_input.fastq.gz,ko_sample2_input,KO,INPUT,pair3;pair4,input,narrow
wt_sample1_chip_prdm.fastq.gz,wt_sample1_chip,WT,PRDM19,pair5,chip,narrow
wt_sample1_chip_h3k4.fastq.gz,wt_sample1_chip,WT,H3K4Me1,pair6,chip,narrow
wt_sample1_input.fastq.gz,wt_sample1_input,WT,INPUT,pair5;pair6,input,narrow
wt_sample2_chip_prdm.fastq.gz,wt_sample2_chip,WT,PRDM19,pair7,chip,narrow
wt_sample2_chip_h3k4.fastq.gz,wt_sample2_chip,WT,H3K4Me1,pair8,chip,narrow
wt_sample2_input.fastq.gz,wt_sample2_input,WT,INPUT,pair7;pair8,input,narrow
```

We have changed file names and descriptions because bcbio does not allow for duplicates. Each chip must be in its own pair ONLY with its associated input file. Since `ko_sample1_chip_prdm` and `ko_sample1_chip_h3k4` are separate chips they get different pairs (pair1 and pair2 respectively) but the input ko_sample1_input matches both pair1 and pair2. Input files can match as many chips as needed and the pair names should be separated by `;`.

The `antibody` column tells bcbio whether to call broad or narrow peaks on the data. Which peak type is the users choice but generally the follow are considered broad antibodies:

    {'h3f3a', 'h3k27me3', 'h3k36me3', 'h3k4me1', 'h3k79me2', 'h3k79me3', 'h3k9me1', 'h3k9me2', 'h4k20me1', 'h3k9me3', 'broad'}

While the following are considered narrow antibodies:

    {'h2afz', 'h3ac', 'h3k27ac', 'h3k4me2', 'h3k4me3', 'h3k9ac', 'narrow'}


If you are not sure which to use, it is best to begin with narrow.


### 2. Generate YAML config file for analysis

```bash
bcbio_nextgen.py -w template metadata/chip-example.yaml metadata/neurons.csv fastq
```

You should see a folder structure:
```
neurons
|---config
|---final
|---work
```

`chip-example/config/neurons.yaml` is the main config file to run the bcbio project. You will
see this file has a copy of the parameters in `chip-example.yaml` for each sample.

### 3. Run the analysis
This will run the analysis on a local machine, using 16 cores.

```bash
cd neurons/work
bcbio_nextgen.py ../config/neurons.yaml -n 16
```

## Parameters

* `peakcaller`: `[macs2]` bcbio just supports MACS2
* `aligner`: supports `bowtie2` and `bwa`. `bwa` will result in a superset of the peaks called by `bowtie2`.
* `keep_duplicates`: do not remove duplicates before peak calling. Defaults to _False_.
* `keep_multimapped`: do not remove multimappers before peak calling. Defaults to _False_.

## Output

### Project directory

```
├── 2023-05-01_neurons
│   ├── ataqv
│   │   ├── index.html -- QC report from ataqv
│   ├── bcbio-nextgen-commands.log -- list of commands run by bcbio
│   ├── bcbio-nextgen.log -- stdout of bcbio log
│   ├── consensus 
│   │   ├── consensus.bed -- consensus peaks from NF fraction
│   │   ├── consensus-counts.tsv -- table of alignments per peak for each sample, calculated by featureCounts
│   ├── data_versions.csv -- versions of data used in the pipeline
│   ├── metadata.csv -- supplied metadata about the samples
│   ├── multiqc
│   │   ├──list_files_final.txt
│   │   ├──multiqc_config.yaml
│   │   ├──multiqc_data
│   │   │  ├──multiqc_bcbio_metrics.txt
│   │   │  ├──multiqc_data_final.json
│   │   │  ├──multiqc_data.json
│   │   │  ├──multiqc_fastqc.txt
│   │   │  ├──multiqc_general_stats.txt
│   │   │  ├──multiqc.log
│   │   │  ├──multiqc_samtools_idxstats.txt
│   │   │  ├──multiqc_samtools_stats.txt
│   │   │  ├──multiqc_sources.txt
│   │   │  ├──seqbuster_isomirs.txt
│   │   │  ├──seqbuster_mirs.txt
│   │   ├── multiqc_report.html -- multiQC report with useful quality control metrics
│   ├──metrics -- folder containing metrics for each sample
│   ├──programs.txt -- versions of programs run in the pipeline
│   ├── project-summary.yaml
```

### Sample directories

```
├── ko_sample1_chip
│   ├── fko_sample1_chip-ready.bam -- all alignments
│   ├── ko_sample1_chip-ready.bam.bai 
│   ├── fko_sample1_chip-ready.bw -- bigwig file of alignments
│   ├── greylist -- info on reads in greylist
│   ├── fastqc -- FASTQC files for the sample and samtools statistics
│   ├── macs2 -- contains peak calls 
```

ready.bam contains only uniquely mapped non-duplicated reads. The stats in the `project/multiqc/multiqc_report.html` include all reads (duplicated, multimappers).

## Downstream analysis

### Quality Control

The **MultiQC** report in the project directory under `multiqc/multiqc_report.html`
has useful quality control information that you can
use to help decide if your ChIP-seq project worked.

It is hard to give specific cutoffs of metrics to use since the kit, the sample
material, the organism, the genome annotations and so on all affect all of the
metrics. We generally look at the samples as a whole for an experiment and see
if any of the samples are outliers in the important metrics. In the **MultiQC**
report, we look at the percentage of reads in the peaks, the mapping percentage,
the [ENCODE library complexity statistics](https://www.encodeproject.org/data-standards/terms/) and the FastQC
metrics to try to spot samples with problems.


#### QC reports

- [MultiQC report](http://atac-userstory.s3-website.us-east-2.amazonaws.com/multiqc_report.html)

### Differential affinity analysis

For doing differential affinity analysis we recommend using
[DiffBind](https://bioconductor.org/packages/release/bioc/html/DiffBind.html). 
The DiffBind tutorials are great for understanding how to go about with your downstream analyses. 

#### Differential affinity reports

- [RMarkdown](http://atac-userstory.s3-website.us-east-2.amazonaws.com/peaks.Rmd)
- [HTML report](http://atac-userstory.s3-website.us-east-2.amazonaws.com/peaks.html)
- [example data](http://atac-userstory.s3-website.us-east-2.amazonaws.com/differential-affinity-example.tar.gz)
