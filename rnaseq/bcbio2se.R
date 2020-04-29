# How to import bcbio project into SummarizedExperiment object
# Run: Rscript bcbio2seq.R [bcbio_project/final/project-dir]

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
