---
title: Tips for bcbio
description: This code helps with fixing jobs which timeout, and other general tips
category: computing
subcategory: bcbio
tags: [bcbio, bash, hpc]
---

## Installing a private bcbio development repository on O2
```bash
wget https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/scripts/bcbio_nextgen_install.py
python bcbio_nextgen_install.py ${HOME}/local/share/bcbio --tooldir=${HOME}/local --nodata
ln -s /n/app/bcbio/biodata/genomes/ ${HOME}/local/share/genomes
mkdir -p ${HOME}/local/share/galaxy
ln -s /n/app/bcbio/biodata/galaxy/tool-data ${HOME}/local/share/galaxy/tool-data
export PATH="${HOME}/local/bin:$PATH"
```

## How to fix potential conda errors during installation
Add the following to your `${HOME}/.condarc`:
```yaml
channels:
  - bioconda
  - defaults
  - conda-forge
safety_checks: disabled
add_pip_as_python_dependency : false
rollback_enabled: false
notify_outdated_conda: false
```

## How to fix jobs bcbio jobs timing out
The O2 cluster can take a really long time to schedule jobs. If you are having problems with bcbio timing out, set your --timeout parameter to something high, like this:
```bash
/n/app/bcbio/dev/anaconda/bin/bcbio_nextgen.py ../config/bcbio_ensembl.yaml -n 72 -t ipython -s slurm -q short -r --tag feany --timeout 6000 -t 0-11:00
```

## How to run a one-node bcbio job (multicore, not multinode)
it just runs a bcbio job on one node of the cluster (no IPython)

```
#!/bin/bash

# https://slurm.schedmd.com/sbatch.html

#SBATCH --partition=priority        # Partition (queue)
#SBATCH --time=3-00:00              # Runtime in D-HH:MM format
#SBATCH --job-name=bcbio            # Job name - any name
#SBATCH -c 10
#SBATCH --mem-per-cpu=10G           # Memory needed per CPU or use mem
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL             # Type of email notification (BEGIN, END, FAIL, ALL)

bcbio_nextgen.py ../config/illumina_rnaseq.yaml -n 10
```

## How to import bcbio project into SummarizedExperiment object

```
library(tidyverse)
library(tximport)
library(SummarizedExperiment)
library(janitor)

project_dir <- args[1]
#project_dir <- "/home/sergey/cluster/2019-08-29_9way_bcbio/final/2019-08-22_project"

metadata <- read_csv(file.path(project_dir, "metadata.csv"))
# clean your metadata if you need
# metadata$batch <- NULL
# metadata$phenotype <- NULL
# colnames(metadata) <- c("sample", "category")
# metadata[2,2] <- "day7"
# metadata[3,2] <- "day7"
# metadata$category <- as.factor(metadata$category)

metrics <- read_tsv(file.path(project_dir, "multiqc", "multiqc_data", "multiqc_bcbio_metrics.txt")) %>% 
    clean_names(case = "snake")

sample_dirs <- file.path(project_dir, "..", metadata$sample)
salmon_files <- file.path(sample_dirs, "salmon", "quant.sf")
names(salmon_files) <- metadata$sample

transcripts2genes_file <- file.path(project_dir, "tx2gene.csv")
transcripts2genes <- read_csv(transcripts2genes_file, col_names = c("ensembl_transcript_id", "ensembl_gene_id"))

txi_salmon <- tximport(salmon_files, type = "salmon", tx2gene = transcripts2genes,
                      countsFromAbundance = "lengthScaledTPM")

raw_counts <- round(data.frame(txi_salmon$counts, check.names = FALSE), 0) %>% as.matrix()

col_data <- metadata %>% column_to_rownames(var = "sample")
col_data$sample <- rownames(col_data)

se_metadata <- list(metrics = metrics,
                 countsFromAbundance = txi_salmon$countsFromAbundance)

vst <- vst(raw_counts)

se <- SummarizedExperiment(assays = list(raw = raw_counts,
                                         tpm = txi_salmon$abundance,
                                         length = txi_salmon$length,
                                         vst = vst),
                                         colData = col_data,
                                         metadata = se_metadata)
saveRDS(se, "bcbio.se.RDS")
```
