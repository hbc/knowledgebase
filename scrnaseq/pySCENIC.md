A collection of SCENIC related questions (whether using R, Python, VSN
or other!) is available
[here](https://github.com/aertslab/SCENIC/discussions).

# Installation

To install `pySCENIC`, see instructions
[here](https://pyscenic.readthedocs.io/en/latest/installation.html).

## Basic installation to run pySCENIC from the CLI

``` r
conda create -y -n pyscenic2 python=3.8.12 
  ## specific python version chosen to match default python loaded by miniconda
conda activate pyscenic2 
  ## Reminder: use `*source* activate` instead of `conda activate` in slurm scripts

conda install -y numpy
conda install -y -c anaconda cytoolz

pip3 install pyscenic
```

## Add-ons for exploration of CLI output in Python

As detailed in “Read The Docs”, there are multiple ways to execute a
`pySCENIC` analysis. This tutorial is made for use of the command line
interface (CLI) and **assumes that all the necessary single-cell RNA-seq
data pre-processing steps were run in R using `Seurat`** (see section
below for data preparation script).

The loom file exported using `SCopeLoomR` seems incompatible with the
latest `loompy` version (installed by default when installing
`pySCENIC`), which prevents the writing of a `pySCENIC` loom output
file. I thus also created a separate environment with an earlier version
of `scanpy` (and default older `loompy`), from where I run the final
steps of the `pySCENIC` analysis.

``` r
conda create -y -n pyloom python=3.7
conda activate pyloom

conda install -y numpy
conda install -y -c anaconda cytoolz

pip install pyscenic
pip install matplotlib 
pip install seaborn
pip install scanpy==1.4.4.post1
  ## iif compatibility issues with loom files

# OPTIONAL:
#conda install -c conda-forge multicore-tsne
  ## only if running dimensional reduction with pySCENIC (I prefer reverting to R/Seurat for this late analysis stage)
```

------------------------------------------------------------------------

# Required databases

Before you run `pySCENIC`, you will further need the following
databases/datasets:

-   A **genome ranking database**, stored in the `feather` format:
    download from
    [cisTarget](https://resources.aertslab.org/cistarget/), or see how
    to [create your
    own](https://github.com/aertslab/create_cisTarget_databases) (I
    haven’t tried that option)
-   A **motif annotation database**, available from “Read The Docs”
    [here](https://pyscenic.readthedocs.io/en/latest/installation.html#auxiliary-datasets)
-   A **list of transcription factors (TFs)** for which you want to find
    regulons (user defined). See literature and/or Human Atlas Protein
    for an exhaustive list. I compiled one from different sources for
    human: see report
    [here](https://github.com/hbc/hbc_medoff_rnaseq_adamts14_human_lung_hbc04302/blob/master/Upstream-TF_siADAMTS14/scripts/medoff_FA_human_siADAMST14_ovlp_gene-annotation.Rmd)
    for the sources used in the compilation of this list, and repo
    [here](https://github.com/hbc/hbc_medoff_rnaseq_adamts14_human_lung_hbc04302/tree/master/Upstream-TF_siADAMTS14/data)
    for the list itself (save as `.tsv` file without header for use with
    `pySCENIC`).

``` r
# Feather databases (recommend submitting as slurm script, as these are large datasets)
wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather
wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather

# Motif annotation database
wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
```

**Important**: Check your `feather` databases were successfully uploaded
using `sha256sum -c sha256sum.txt` (file provided on
[cisTarget](https://resources.aertslab.org/cistarget/) webpage)

------------------------------------------------------------------------

# Single cell data preparation (using R)

To run `pySCENIC` from the command line interface, all you need is a
loom file of your single-cell experiment. This can be prepared from
`scanpy` (the `pySCENIC` tutorial includes an example of pre-processing
with `scanpy`
[here](http://htmlpreview.github.io/?https://github.com/aertslab/SCENICprotocol/blob/master/notebooks/PBMC10k_SCENIC-protocol-CLI.html))
or with `Seurat`.

I used `Seurat` and exported my processed (QC’d, normalized, integrated,
and clustered) object as a loom file from R, as described in the
[`SCopeLoomR`
vignette](https://htmlpreview.github.io/?https://github.com/aertslab/SCopeLoomR/blob/master/vignettes/SCopeLoomR_tutorial.nb.html).
This requires the `SCopeLoomR` package, installed using
`devtools::install_github("aertslab/SCopeLoomR")`.

**The steps detailed in this section are wrapped into one single R
script - see example from the [Aikawa
consult](https://github.com/hbc/aikawa_scRNA-scATAC-seq_human_macrophage_hiv-vesicles_hbc04403/blob/main/R_O2/scripts/5C_aikawa_pyscenic_seurat_pre-processing.R).**

## Seurat object upload

First, upload your `Seurat` object and extract your raw (or normalized)
count matrix from it.

``` r
library(SCopeLoomR)

### Load minimal Seurat object with RNA assay
seurat_rna <- readRDS("../data/3A_seurat_RNA_2022-04-01_combined.rds")

# Extract raw count matrix [can also use normalized]
dgem <- seurat_rna@assays$RNA@counts
dim(dgem)
head(colnames(dgem)) # columns => cell IDs

# Extract cell-level metadata
cell.info <- seurat_rna@meta.data
```

## Count matrix filtering

**This pre-filtering step seems critical.** All the jobs I tried to run
stalled without error messages until I ran this step. Here, I restricted
the analysis to genes detected in at least 1 per thousand cells, and
further excluded mitochondrial and ribosomal genes.

``` r
# Compute logical vector for every gene reflecting whether
# it has more than zero count per cell & and is detected in at 1 per thousand of cells in the global dataset
c <- floor(dim(dgem)[2]*.001); c  # c = 50 cells
nonzero <- dgem > 0L
keep_genes <- rowSums(as.matrix(nonzero)) >= c
filtered_dgem <- dgem[keep_genes, ]

# Further exclude mitochondrial genes
idx <- grep("^MT-", rownames(filtered_dgem))
rownames(filtered_dgem)[idx]
filtered_dgem <- filtered_dgem[-idx, ]
dim(filtered_dgem) # 20320 genes x 50931 cells

# Further exclude ribosomal genes [NA for this dataset]
idx <- grep("^RP.-", rownames(filtered_dgem))
rownames(filtered_dgem)[idx]
#filtered_dgem <- filtered_dgem[-idx, ]
#dim(filtered_dgem)
```

## Create loom file

Once the count matrix is ready, you can generate the loom file. It is
possible to add as much information as you want to the loom file
(embeddings, cluster information, etc.), though I would recommend to
keep the file as light as possible and do all the post-processing steps
in R.

``` r
### OPTIONAL (necessary iif using pySCENIC for dimensionality analyses)

# Extract default embedding (e.g. UMAP or PCA coordinates)
default.umap <- Embeddings(seurat_rna, reduction = "umap")
default.umap.name <- "UMAP"


### Create loom file

file.name <- "3A_seurat_RNA_2022-04-01_combined.loom"
project.title <- "Aikawa multiomics complete dataset [2022-04-01]"
build_loom(
  file.name = file.name,
  dgem = dgem,
  title = project.title,
  genome = "human"#,
  #default.embedding = default.umap,
  #default.embedding.name = default.umap.name
)


### OPTIONAL

# Add further information to loom file (e.g. other embeddings)
loom <- open_loom(file.path = file.name, mode = "r+")

other.umap <- Embeddings(seurat_rna, reduction = "umap.rna")
add_embedding(loom = loom, embedding = other.umap, name = "UMAPxRNA")

pca.reduc <- Embeddings(seurat_rna, reduction = "pca")
add_embedding(loom = loom, embedding = other.umap, name = "PCA")

# Add Seurat cluster definition
add_seurat_clustering(loom = loom, seurat = seurat_rna,
                      seurat.clustering.prefix = "wsnn_res.", 
                      default.clustering.resolution = "res.0.15")
# can also add cluster markers and other cluster-level annotations (see vignette)

finalize(loom = loom)  # close loom file after add-ons
```

------------------------------------------------------------------------

# Run CLI implementation of pySCENIC

**Note**: It is also possible to run `pySCENIC` from Python (as
described
[here](https://pyscenic.readthedocs.io/en/latest/tutorial.html#zeisel-et-al-dataset)).

For easiness of use in the O2 HPC environment, I used the command line
implementation (i.e. extracting the minimal CLI commands from this
[vignette](http://htmlpreview.github.io/?https://github.com/aertslab/SCENICprotocol/blob/master/notebooks/PBMC10k_SCENIC-protocol-CLI.html),
which also explains how to pre-process the raw single-cell data with
`scanpy`).

**Steps 1 and 2 are wrapped into one single slurm script - see example
from the [Aikawa
consult](https://github.com/hbc/aikawa_scRNA-scATAC-seq_human_macrophage_hiv-vesicles_hbc04403/blob/main/slurm_scripts/5_aikawa_pySCENIC_part1.sbatch).**
Step 3 is run within a second slurm script, as it requires fewer
resources and a different conda environment - see example from the
[Aikawa
consult](https://github.com/hbc/aikawa_scRNA-scATAC-seq_human_macrophage_hiv-vesicles_hbc04403/blob/main/slurm_scripts/5_aikawa_pySCENIC_part2.sbatch).

## Step 1: Find co-expression networks

### Using `pyscenic grn`

For small datasets (= 2k cells x 10k genes), I used the following
script, which takes about 2 hours to run:

``` r
#!/bin/bash

#SBATCH -p priority             # partition name
#SBATCH -c 16                   # Request N cores (must match future plan requested in R script)
#SBATCH -t 0-6:00               # hours:minutes runlimit after which job will be killed
#SBATCH --mem 96G               # total memory
#SBATCH --job-name pyscenic     # Job name
#SBATCH -o %j.out               # File to which standard out will be written
#SBATCH -e %j.err               # File to which standard err will be written


# Move to desired working directory
cd /n/data1/cores/bcbio/PIs/masanori_aikawa/
cd aikawa_scRNA-scATAC-seq_human_macrophage_hiv-vesicles_hbc04403/pySCENIC

# Set up Python environment for pySCENIC
module load conda2
source activate pyscenic

# Create results directory
mkdir -p outs_2022-05-02/

# Run pySCENIC get regulatory network from command line interface
pyscenic grn \
        --num_workers 16 \
        -o outs_2022-05-02/3A_seurat_RNA_pySCENIC_GRN_adjacencies.csv \
        ../R_O2/data/3A_seurat_RNA_2022-04-01_combined.loom \
        scenic_db/human_transcription-factors_curated_list.tsv

echo "regulatory networks defined"
```

### Using `arboreto_with_multiprocessing.py`

With larger datasets (>10k cells x 20k genes), I had trouble with `dask`
(workers timing out in the middle of the job and leading to errors,
similar to what users describe
[here](https://github.com/aertslab/pySCENIC/issues/54)).

Therefore, I used the following approach instead, as recommended in the
`pySCENIC`
[FAQ](https://pyscenic.readthedocs.io/en/latest/faq.html#i-am-having-problems-with-dask).

``` r
#!/bin/bash

#SBATCH -p short                # partition name
#SBATCH -c 18                   # Request N cores [>12 cores recommended for pyscenic2 steps]
#SBATCH -t 0-8:00               # hours:minutes time after which job will be killed [16h needed for 50k cells x 20k cells with 20 cores]
#SBATCH --mem 120G              # total amount of memory requested [>100G recommended]
#SBATCH --job-name henderson    # Job name
#SBATCH -o %j.out               # File to which standard out will be written
#SBATCH -e %j.err               # File to which standard err will be written


# Move to desired working directory
cd /home/aj186/10X_oJIA_v3_re-analysis/SF_Teff+Treg

# Set up Python environment
module load miniconda3/4.10.3
  ## This module is necessary to use *source* activate pyscenic2
source activate pyscenic2

# Double-check what Python version is running
python --version
echo "PATH ="; echo $PATH
echo "PYTHONPATH ="; echo $PYTHONPATH


# Create results directory
res_dir=outs_$(date +"%Y-%m-%d")
mkdir -p $res_dir
date

# Run pySCENIC get regulatory network from command line interface
# This circumvents the use of dask, which seems to cause trouble...
arboreto_with_multiprocessing.py \
        /home/aj186/10X_oJIA_v3_re-analysis/SF_Teff+Treg/data/3_integrated_Seurat_SF_Teff+Treg_filtered.loom \
        scenic_db/human_transcription-factors_curated_list.tsv \
        --num_workers 18 \
        -o $res_dir/3A_integrated_SF_Teff+Treg_pySCENIC_GRN_adjacencies.csv \
        --method grnboost2 \
        --seed 787878

date
echo "regulatory networks inferred"
```

## Step 2: Define modules of TF-regulons

Directly continuing from the same slurm script as above, I run
`pyscenic ctx`. This step cleans up the identified co-expression
networks by drawing information from the feather databases and motif
annotation files to identify target genes with significantly enriched
motifs.

``` r
# Clean up modules using feather ranking databases
pyscenic ctx $res_dir/3A_integrated_SF_Teff+Treg_pySCENIC_GRN_adjacencies.csv \
    scenic_db/*.mc9nr.feather \
    --annotations_fname scenic_db/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
    --expression_mtx_fname /home/aj186/10X_oJIA_v3_re-analysis/SF_Teff+Treg/data/3_integrated_Seurat_SF_Teff+Treg_filtered.loom \
    --output $res_dir/3A_integrated_SF_Teff+Treg_pySCENIC_CTX_regulons.csv \
    --mask_dropouts \
    --num_workers 18

date
echo "regulon modules defined"
sleep 5
```

## Step 3: Calculate regulon modules score for each cell

Finally, I run `pyscenic aucell` to calculate a regulon module score for
each cell in the dataset. This step is faster and much less
resource-intensive than the steps above. This is also the step where I
need to switch to my pyloom environment, otherwise the loom output files
do not get written.

``` r
#!/bin/bash

#SBATCH -p short                # partition name
#SBATCH -c 12                   # Request N cores (must match future plan requested in R script)
#SBATCH -t 0-8:00               # hours:minutes runlimit after which job will be killed
#SBATCH --mem 96G               # total amount of memory requested
#SBATCH --job-name henderson    # Job name
#SBATCH -o %j.out               # File to which standard out will be written
#SBATCH -e %j.err               # File to which standard err will be written

# Move to desired working directory
cd /home/aj186/10X_oJIA_v3_re-analysis/SF_Teff+Treg

# Set up Python environment
module load miniconda3/4.10.3
source activate pyloom
  ## Reminder: This environment uses older versions of Python/packages (esp. loom)
  ##           This is necessary for compatibility with loom file v2 generated from R

# Double-check what Python version is running
python --version
echo "PATH ="; echo $PATH
echo "PYTHONPATH ="; echo $PYTHONPATH


# Go to results directory from pyscenic run (part 1)
res_dir=outs_2022-05-16
date

# Calculate cell module scores (from filtered loom, to extract acurate regulon incidence matrix later on)
pyscenic aucell \
    /home/aj186/10X_oJIA_v3_re-analysis/SF_Teff+Treg/data/3_integrated_Seurat_SF_Teff+Treg_filtered.loom \
    $res_dir/3A_integrated_SF_Teff+Treg_pySCENIC_CTX_regulons.csv \
    --output $res_dir/3A_integrated_SF_Teff+Treg_2022-05-16_pySCENIC_filtered.loom \
    --num_workers 12

echo "cell scoring completed"

# Add dimensionality reduction (derived from pySCENIC AUCell matrix)
   ## I always skipped this step, as it can also be performed in R (with proper sample integration)
#python add_visualization.py \
#    --loom_input $res_dir/3A_seurat_RNA_2022-05-12_pySCENIC.loom \
#    --loom_output $res_dir/3A_seurat_RNA_2022-05-12_pySCENIC_viz.loom \
#    --num_workers 20

#echo "dimensionality reduction added"
```

### Note on final step (dimensionality reduction with pySCENIC)

To enable visualization in `SCope` and other downstream analyses, it can
be useful to have a dimensionality reduction based on the regulon scores
for each cell. This is not automatically added when running `AUCell`
from the CLI, but can be added using the following python script
[`add_visualization.py`](https://github.com/vib-singlecell-nf/vsn-pipelines/tree/58137baa31e580e82e5e2add59e2b536d9754bd0/src/scenic/bin).

**Although “Read The Docs” suggests this script is automatically added
when installing `pySCENIC`, this was not the case for me, and I add to
grab the Python script from the GitHub repo linked above before I could
run it.** (This is another reason why I prefer to revert to R and
`Seurat` for these downstream analysis steps.)

------------------------------------------------------------------------

# Integration of pySCENIC output with Seurat

From there, I use the following R script to extract relevant outputs
from the `pySCENIC` loom file and integrate them with my initial
`Seurat` object. **The steps detailed in this section are wrapped into
one single R script - see example from the [Aikawa
consult](https://github.com/hbc/aikawa_scRNA-scATAC-seq_human_macrophage_hiv-vesicles_hbc04403/blob/main/R_O2/scripts/5C_aikawa_pyscenic_seurat_post-processing.R).**

``` r
# Set up R environment
library(Seurat)
library(Signac)
library(tidyverse)
library(data.table)
library(SCopeLoomR)


### Define study-specific parameters

out_dir <- "outs_2022-05-16"
out_files <- list.files(out_dir, full.names = TRUE)

seurat_path <- "/home/aj186/10X_oJIA_v3_re-analysis/SF_Teff+Treg/data/3_integrated_Seurat_SF_Teff+Treg.rds"
analysis <- "3_integrated_SF_Teff+Treg"


### Load pySCENIC output loom file (regulons)

# from pySCENIC working directory:

path2loom <- out_files[grep("_pySCENIC_filtered.loom", out_files)]
  # using *filtered* loom to extract regulons
loom <- open_loom(path2loom, mode = "r")

# Read pySCENIC information from loom file
# WARNING: column attribute names may be different depending on the SCopeLoomR/pySCENIC version used!

# => regulons
regulonsMat <- get_regulons(loom,  # as incidence matrix
                            column.attr.name = "Regulons")
dim(regulonsMat)
write.csv(regulonsMat,
          file = paste0(out_dir, "/", analysis, "_pySCENIC_regulons_incidence_matrix.csv"))

regulons <- SCENIC::regulonsToGeneLists(regulonsMat)  # can convert to gene list for processing with R
head(regulons)

close_loom(loom)


### Add desired data to original Seurat object (RNA only)

seurat_rna <- readRDS(seurat_path)

# Extract underlying matrix from regulonsAUC object
AUCmat <- AUCell::getAUC(regulonsAUC)
AUCmat[1:10, 1:10]
rownames(AUCmat)

# Add as assay to Seurat object
seurat_rna[['pyscenicAUC']] <- CreateAssayObject(data = AUCmat)
saveRDS(seurat_rna, file = paste0(out_dir, "/", analysis, "_pyscenic2seurat.rds"))
```

After running that script, any `Seurat` tool can be used for
dimensionality reduction and re-clustering by regulon score, or for
calculation of differential activity scores. I do not detail these steps
here, but see for example this R markdown report from the [Aikawa
consult](https://github.com/hbc/aikawa_scRNA-scATAC-seq_human_macrophage_hiv-vesicles_hbc04403/blob/main/R_O2/scripts/5_aikawa_cluster-GRN_pySCENIC_human_macrophage_multiomics.Rmd).

The `pySCENIC` documentation also includes a vignette for
post-processing and interpretation using `scanpy`
[here](http://htmlpreview.github.io/?https://github.com/aertslab/SCENICprotocol/blob/master/notebooks/PBMC10k_downstream-analysis.html).
