# ChIP-seq

bcbio evaluates data quality using FASTQC and filters and trims reads as necessary. Reads are aligned to the genome and only uniquely aligning reads are retained. Peaks are called with MACS2. A consensus peak file is generated using bedops and this consensus peak list is used to generate a count matrix using featureCounts.

## Description of example dataset

We will be working using ChIP-seq data from a recent publication in Neuron by Baizabal et al. (2018) [[1]](https://doi.org/10.1016/j.neuron.2018.04.033). 
Baizabal et al. sought to understand how chromatin-modifying enzymes function in neural stem cells to establish the epigenetic landscape that determines cell type and stage-specific gene expression. Chromatin-modifying enzymes are transcriptional regulators that control gene expression through covalent modification of DNA or histones. 

This specific dataset focuses on the transcriptional regulator **PRDM16, which is a chromatin-modifying enzyme** that belongs to the larger PRDM (Positive Regulatory Domain) protein family, that is structurally defined by the **presence of a conserved N-terminal histone methyltransferase PR domain** ([Hohenauer and Moore, 2012](https://journals.biologists.com/dev/article/139/13/2267/45169/The-Prdm-family-expanding-roles-in-stem-cells-and)). The authors generated CHiP-seq data for PRDM16. Our dataset consists of two WT samples and two KO samples.


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

#### 1.2 Copy the template YAML file describing the CHiP-seq analysis


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


```
File,description,genotype,enzyme,batch,phenotype,antibody
ko_sample1_chip.fastq.gz,ko_sample1_chip,KO,PRDM19,pair1,chip,narrow
ko_sample1_input.fastq.gz,ko_sample1_input,KO,PRDM19,pair1,input,narrow
ko_sample2_chip.fastq.gz,ko_sample2_chip,KO,PRDM19,pair2,chip,narrow
ko_sample2_input.fastq.gz,ko_sample2_input,KO,PRDM19,pair2,input,narrow
wt_sample1_chip.fastq.gz,wt_sample1_chip,WT,PRDM19,pair3,chip,narrow
wt_sample1_input.fastq.gz,wt_sample1_input,WT,PRDM19,pair3,input,narrow
wt_sample2_chip.fastq.gz,wt_sample2_chip,WT,PRDM19,pair4,chip,narrow
wt_sample2_input.fastq.gz,wt_sample2_input,WT,PRDM19,pair4,input,narrow
```

The necessary columns here are: `File`, `description`, `batch`, `phenotype` and `antibody`. 
For ChIP-seq, bcbio requires `batch`,`phenotype`, and `antibody` are unique to ChIP-seq.

`batch` matches your input samples with their respective chips and the `phenotype` column tells bcbio if a sample is an input or chip.

Here we have one input for every chip. For example ko_sample1 has ko_sample1_chip and ko_sample1_input. These are pair1. However, sometimes the same input is used for multiple chips. Here is the same file but assuming that we also ran a `h3k4me1` chip on all samples


```

```


However, please note that the `antibody` column should be added with caution.

- Valid broad antibodies are: 

    {'h3f3a', 'h3k27me3', 'h3k36me3', 'h3k4me1', 'h3k79me2', 'h3k79me3', 'h3k9me1', 'h3k9me2', 'h4k20me1', 'h3k9me3', 'broad'}

- Valid narrow antibodies are: 

    {'h2afz', 'h3ac', 'h3k27ac', 'h3k4me2', 'h3k4me3', 'h3k9ac', 'narrow'}


If you know your antibody should be called with narrow or broad peaks, supply 'narrow' or 'broad' as the antibody.
```
Bcbio will call narrow peaks if you have a antibody column, but do not have a vaild antibody within that list.

By default, you will get *.narrowPeak files if you do not have a antibody column.
```

### 2. Generate YAML config file for analysis
```bash
bcbio_nextgen.py -w template metadata/atac-example.yaml metadata/hindbrain_forebrain.csv fastq
```

In the result you should see a folder structure:
```
hindbrain_forebrain
|---config
|---final
|---work
```

`hindbrain_forebrain/config/hindbrain_forebrain.yaml` is the main config file to run the bcbio project. You will
see this file has a copy of the parameters in `atac-example.yaml` for each sample.

### 3. Run the analysis
This will run the analysis on a local machine, using 16 cores.
```bash
cd hindbrain_forebrain/work
bcbio_nextgen.py ../config/hindbrain_forebrain.yaml -n 16
```

## Parameters
* `peakcaller`: `[macs2]` bcbio just supports MACS2
* `aligner`: supports `bowtie2` and `bwa`. `bwa` will result in a superset of the peaks called by `bowtie2`.
* `chip_method`: set to `atac` to run the ATAC-seq pipeline
* `keep_duplicates`: do not remove duplicates before peak calling. Defaults to _False_.
* `keep_multimapped`: do not remove multimappers before peak calling. Defaults to _False_.

## Output

### Project directory

```
├── 2020-05-01_hindbrain_forebrain
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
│   │   ├── multiqc_report.html -- multiQC report with useful quality control metrics
│   ├── programs.txt -- versions of programs run in the pipeline
```

### Sample directories

```
├── forebrain_rep1
│   ├── forebrain_rep1-DN.bam -- dinucleosome alignments
│   ├── forebrain_rep1-full.bam -- all fraction alignments
│   ├── forebrain_rep1-MN.bam -- mononucleosome alignments
│   ├── forebrain_rep1-NF.bam -- nucleosome-free alignments
│   ├── forebrain_rep1-ready.bam -- identifical to -full
│   ├── forebrain_rep1-ready.bam.bai 
│   ├── forebrain_rep1-ready.bw -- bigwig file of full alignments
│   ├── forebrain_rep1-TN.bam -- trinucleosome alignments
│   ├── macs2 -- contains peak calls for each fraction, including the full peak calls
```

read.bam contains only uniquely mapped non-duplicated reads, see [bam cleaning function](https://github.com/bcbio/bcbio-nextgen/blob/master/bcbio/chipseq/__init__.py#L18). The stats in the `project/multiqc/multiqc_report.html` include all reads (duplicated, multimappers).

## Downstream analysis

### Quality Control
The **MultiQC** report in the project directory under `multiqc/multiqc_report.html`
and at the **ataqv** report in the project directory under
`ataqv/ataqv_report.html` have useful quality control information that you can
use to help decide if your ATAC-seq project worked.

It is hard to give specific cutoffs of metrics to use since the kit, the sample
material, the organism, the genome annotations and so on all affect all of the
metrics. We generally look at the samples as a whole for an experiment and see
if any of the samples are outliers in the important metrics. In the **MultiQC**
report, we look at the percentage of reads in the peaks, the mapping percentage,
the 
[ENCODE library complexity statistics](https://www.encodeproject.org/data-standards/terms/) and the FastQC
metrics to try to spot samples with problems.

In the **ataqv** report, we look at the HQAA fragment length distribution plot.
Ideally, this plot should show a periodic uptick every 200 bases, which
corresponds to the different nucleosome fractions. The samples should be
enriched for < 100 which is the nucleosome free fraction, 200 for the
mononucleosome fraction, 400 for the dinucleosome fraction and 600 for the
trinucleosome fraction. Often you will not see this behavior though even in
libraries that were successful. But if some of your samples have this and others
do not, that is something to be concerned about.

You should see an enrichment around the transcription start sites, if you are
missing that then your experiment likely failed. The **peaks** table in the
**tables** tab in the **ataqv** report has a measurement of the high quality
autosomal alignments overlapping peaks, **ataqv** calculates this metric using
all of the peaks, not just the peaks from the nucleosome-free fraction, so this
is useful to look at as well. See the [ataqv github
repository](https://github.com/ParkerLab/ataqv/issues/13) for a discussion of
the ranges of values you can expect to see for metrics in the **ataqv** report
along with other values to look at that might be informative. The Parker lab
reprocessed samples from many publications with **ataqv** and posted the reports
[here](https://theparkerlab.med.umich.edu/data/porchard/ataqv-public-survey/)
which is helpful to browse through to get an idea of what ranges of values you
can expect. As you can see, they can be all over the place.

#### hindbrain vs forebrain QC reports
- [MultiQC report](http://atac-userstory.s3-website.us-east-2.amazonaws.com/multiqc_report.html)
- [ataqv report](http://atac-userstory.s3-website.us-east-2.amazonaws.com)

### Differential affinity analysis
For doing differential affinity analysis we recommend loading the consensus peak
table from the `consensus/consensus.counts` file in the project directory. This
is a table of counts per peak in the nucleosome-free fraction for each sample
that you can use in any standard count-based differential expression tools like
[DESeq2](
https://bioconductor.org/packages/release/bioc/html/DESeq2.html)/[edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)/[limma](https://bioconductor.org/packages/release/bioc/html/limma.html).
Often you will find tutorials using
[DiffBind](https://bioconductor.org/packages/release/bioc/html/DiffBind.html)
but that uses these callers under the hood, so you can just call them directly
and skip an intermediate step if you want. Either way works. The DiffBind tutorials
are great for understanding how to go about with your downstream analyses.

#### hindbrain vs forebrain differential affinity reports
- [RMarkdown](http://atac-userstory.s3-website.us-east-2.amazonaws.com/peaks.Rmd)
- [HTML report](http://atac-userstory.s3-website.us-east-2.amazonaws.com/peaks.html)
- [example data](http://atac-userstory.s3-website.us-east-2.amazonaws.com/differential-affinity-example.tar.gz)
