# Leafcutter Splice Variant Analysis

Leafcutter documentation is found [HERE](http://davidaknowles.github.io/leafcutter/))

### About Leafcutter

Leafcutter quantifies RNA splicing variation using short-read RNA-seq data. The core idea is to leverage spliced reads (reads that span an intron) to quantify (differential) intron usage across samples. 
LeafCutter uses short-read RNA-seq data to detect intron excision events at base-pair precision by analyzing mapped split reads. LeafCutter focuses on alternative splicing events, including skipped exons, 5′ and 3′ alternative splice-site usage, and additional complex events that can be summarized by differences in intron excision. LeafCutter’s intron-centric view of splicing is based on the observation that mRNA splicing occurs predominantly through the step-wise removal of introns from nascent pre-mRNA. (Unlike isoform-quantification methods such as Cufflinks2, LeafCutter does not measure alternative transcription start sites and alternative polyadenylation directly, as they are not generally captured by intron excision events.) The major advantage of this representation is that LeafCutter does not require read assembly or inference of isoforms supported by ambiguous reads, both of which are computationally and statistically difficult.

### Installing Leafcutter

#### Prerequisites

samtools should be available on your PATH
regtools should be available on your PATH - regtools is easy to [install](https://regtools.readthedocs.io/en/latest/#installation)
Python 2.7 (earlier versions may be OK)
R (version 3.6.0, earlier versions may be OK)

**Leafcutter requires the R package rstan which is incredibly difficult to install. Here are instructions from HMS (note that installation takes approx 1 hour:**

1. Start an interactive session with a lot of memory
```bash
srun --pty -p interactive -t 0-3 --mem 10G bash
```

2. Purge any loaded modules and re-install. This installation was for R 4.1.1

```bash 
module purge 
module load gcc/6.2.0 R/4.1.1 
```
3. Make a home library and export it. If you already have a library for R you can skip the first command

```bash
mkdir -p ~/R-4.1.1-rstan/library 
echo 'R_LIBS_USER="~/R-4.1.1-rstan/library"' > $HOME/.Renviron 
export R_LIBS_USER="~/R-4.1.1-rstan" 
```

4. Start R and make sure that you remove any prior attempts to install rstan

```R
R 
> remove.packages("rstan") 
> if (file.exists(".RData")) file.remove(".RData") 
> quit() # do not save workspace image
```

5. Restart R and install R stan. Then test the installation.

```R
R 
> Sys.setenv(DOWNLOAD_STATIC_LIBV8 = 1) 
> install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE) # this part will take about 45 mins
> example(stan_model, package = "rstan", run.dontrun = TRUE) # test the installation
```

#### Installing leafcutter itself

Once you have done all the pre-reqs leafcutter itself is incredibly easy to install in R

```R
devtools::install_github("davidaknowles/leafcutter/leafcutter")
```

## Running Leafcutter

I ran leafcutter on an interactive session on O2 with a lot of memory

```bash
srun --pty -p interactive -t 0-5 --mem 10G bash
```

### Step 1

Step 1 converts bam to junction files. Your chromosome names must be chr1...chrN plus chrX and chrY. regtools will skip chromosomes named otherwise.

I ran this script iteractively

```bash
#!/bin/sh

module load gcc/9.2.0
module load samtools/1.14

for bamfile in `ls /n/data1/cores/bcbio/PIs/peter_sicinski/sicinski_inhibition_RNAseq_human_hbc04676/final/*/*ready.bam`; do
    echo Converting $bamfile to $bamfile.junc
samtools index $bamfile

/n/data1/cores/bcbio/PIs/peter_sicinski/sicinski_inhibition_RNAseq_human_hbc04676/leafcutter/regtools/build/regtools junctions extract -a 8 -m 50 -s RF  -M 500000 $bamfile -o $bamfile.junc 
    echo $bamfile.junc >> test_juncfiles.txt
done
```

test_juncfiles.txt just has the locations of the junction files. I had three different experiments so split this file into three files before proceeding to step 2.

### Step 2 

Step 2 is intron clustering. I ran this three times, once for each experiment (TD47, NCIH_PD, and NCIH_SRPIN)

```bash
python leafcutter/clustering/leafcutter_cluster_regtools.py -j T47D_juncfiles.txt -m 50 -o TD47  -l 500000
python leafcutter/clustering/leafcutter_cluster_regtools.py -j NCIH_PD_juncfiles.txt -m 50 -o NCIH_PD  -l 500000
python leafcutter/clustering/leafcutter_cluster_regtools.py -j NCIH_SRPIN_juncfiles.txt -m 50 -o NCIH_SRPIN  -l 500000
```

### Step 3

Step three test for differential intron usage. To run step 3 you need a groups file which is basically a metadata file. Here is the one that I used for my TD47 dataset. These are tab separated files.

-----------------------------------------------------

T47D_K6_D3_overexpression_replicate_3-ready.bam	K6   
T47D_empty_vector_replicate_2-ready.bam	empty    
T47D_empty_vector_replicate_1-ready.bam	empty    
T47D_empty_vector_replicate_3-ready.bam	empty     
T47D_K6_D3_overexpression_replicate_2-ready.bam	K6     
T47D_K6_D3_overexpression_replicate_1-ready.bam	K6      

------------------------------------------------------

The actual analysis was run on O2 in the same interactive session as before.Note that the -i 3 parameter was necessary because I only had 3 samples per group not 4 (the default minimum). The exon file was downloaded from the leafcutter github on March 28, 2023 https://www.dropbox.com/s/osmtksukti9gv7o/gencode.v31.exons.txt.gz?dl=0 


I ran the code 3x, once per data set

```bash
Rscript leafcutter/scripts/leafcutter_ds.R --num_threads 1 NCIH_SRPIN_perind_numers.counts.gz SRPIN_groups.txt -i 3 -e gencode.v31.exons.txt.gz 

Rscript leafcutter/scripts/leafcutter_ds.R --num_threads 1 NCIH_PD_perind_numers.counts.gz PD_groups.txt -i 3 -e gencode.v31.exons.txt.gz 

Rscript leafcutter/scripts/leafcutter_ds.R --num_threads 1 TD47_perind_numers.counts.gz T47D_groups.txt -i 3 -e gencode.v31.exons.txt.gz 
```


### Merging output in R

Leafcutter outputs 2 separate files for each contrast you run. leafcutter_ds_cluster_significance_CONTRASTNAME.txt and leafcutter_ds_effect_sizes_CONTRASTNAME.txt. Below is R code to merge these two files and summarize at the level of intron cluster. 

Note that this script takes 1 parameter (psi_cutoff) which is 0.01 in the below script but can be modified. The script will count the number of introns in the cluster above this value.


```R
#determine parameter value
psi_cutoff = 0.01


#read in p-values and annotate

pval_K6 <- read.table("/path/leafcutter_ds_cluster_significance_PD.txt",sep="\t", header=TRUE) %>%
  left_join(pruned_annotations, 
            by = c("genes" = "external_gene_name"))   ## Pruned annotations is a biomart annotations dataframe


#only keep good runs (status=success)

pval_PD <- subset(pval_PD, pval_PD$status=="Success")
pval_PD$N_introns <- (pval_PD$df +1)
effect_PD <- read.table("/path/leafcutter_ds_effect_sizes_PD.txt",sep="\t", header=TRUE)
effect_PD$cluster = paste0(sapply( str_split(effect_PD$intron, ":"),"[", 1 ),":",sapply(str_split(effect_PD$intron, ":"),"[", 4 )) #add cluster name

effect_PD_sum <- effect_PD %>%
 group_by(cluster) %>% summarize(count_deltapsi_above_cutoff = sum(abs(deltapsi) >psi_cutoff)) #determine how many introns are above deltapsi cutoff
effect_PD_max <- effect_PD %>%
  group_by(cluster) %>% summarise_each(funs(deltapsi[which.max(abs(deltapsi))])) # determine maximum value of the deltapsi

effect_PD_summary = effect_PD_sum %>% left_join(effect_PD_max, by=c("cluster")) #merge to have all info

LC_full_PD_info <- pval_PD %>% left_join(dplyr::select(effect_PD_summary, "cluster","count_deltapsi_above_cutoff","deltapsi"), by=c("cluster"="cluster")) # merge again with p-values

colnames(LC_full_PD_info)[14] <- "max_deltapsi"

```
