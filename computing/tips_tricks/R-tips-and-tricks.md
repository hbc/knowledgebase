---
title: R tips
description: This code helps with regular data improving efficiency
category: computing
subcategory: tips_tricks
tags: [R, visualization]
---

# Import/Export of files
Stop using write.csv, write.table and use the [rio](https://cran.r-project.org/web/packages/rio/index.html) library instead. All rio needs is the file extension to figure out what file type you're dealing with. Easy import and export to Excel files for clients.

# Parsing in R using Tidyverse
This is a link to a nice tutorial from Ista Zahn from IQSS using stringr and tidyverse for parsing files in R. It is from the Computefest 2017 workshop:
http://tutorials-live.iq.harvard.edu:8000/user/zwD2ioESyGbS/notebooks/workshops/R/RProgramming/Rprogramming.ipynb

# Better clean default ggplot
install cowplot (https://cran.r-project.org/web/packages/cowplot/index.html)
```r
library(cowplot)
```

# Nice looking log scales
Example for x-axis
```r
library(scales)
   p + scale_x_log10(
         breaks = scales::trans_breaks("log10", function(x) 10^x),
         labels = scales::trans_format("log10", scales::math_format(10^.x))) +
       annotation_logticks(sides='b')
```

# Read a bunch of files into one dataframe
```r
library(purrr)
library(readr)
library(dplyr)
library(tidyr)
read_files = function(files) {
  data_frame(filename = files) %>%
    mutate(contents = map(filename, ~ read_tsv(.))) %>%
    unnest()
}
```

# remove a layer from a ggplot2 object with ggedit
```
plotGeneSaturation(bcb, interestingGroups=NULL) +
  ggrepel::geom_text_repel(aes(label=description, color=NULL))
p %>%
  ggedit::remove_geom('point', 1) +
  geom_point(aes(color=NULL))
```

# [Link to information about count normalization methods](https://github.com/hbc/knowledgebase/wiki/Count-normalization-methods)
The images currently break, but I will update when the course materials are in a more permanent state.

# .Rprofile usefulness
```R
## don't ask for CRAN repository
options("repos" = c(CRAN = "http://cran.rstudio.com/"))
## for the love of god don't open up tcl/tk ever
options(menu.graphics=FALSE)
## set seed for reproducibility
set.seed(123456)
## don't print out more than 100 lines at once
options(max.print=100)
## helps with debugging Bioconductor/S4 code
options(showErrorCalls = TRUE, showWarnCalls = TRUE)

## Create a new invisible environment for all the functions to go in
## so it doesn't clutter your workspace.
.env <- new.env()

## ht==headtail, i.e., show the first and last 10 items of an object
.env$ht <- function(d, n=10) rbind(head(d, n), tail(d, n))

## copy from clipboard
.env$pbcopy = function(x) {
  capture.output(x, file=pipe("pbcopy"))
}

## update your local bcbioRNASeq and bcbioSingleCell installations
.env$update_bcbio = function(x) {
    devtools::install_github("steinbaugh/basejump")
    devtools::install_github("hbc/bcbioBase")
    devtools::install_github("hbc/bcbioRNASeq")
    devtools::install_github("hbc/bcbioSingleCell")
}

attach(.env)
```

# Make density plot without underline
```R
ggplot(colData(sce) %>%
       as.data.frame(), aes(log10GenesPerUMI)) +
    stat_density(geom="line") +
    facet_wrap(~period + intervention)
```

# Archive a file to Dropbox with a link to it
```R
```{r results='asis'}
dropbox_dir = "HSPH/eggan/hbc02067"
archive_data_with_link = function(data, filename, description, dropbox_dir) {
    readr::write_csv(data, filename)
    links = bcbioBase::copyToDropbox(filename, dropbox_dir)
    link = gsub("dl=0", "dl=1", links[[1]]$url)
    basejump::markdownLink(filename, link, paste0(" ", description))
}
archive_data_with_link(als, "dexseq-all.csv", "All DEXSeq results", dropbox_dir)
archive_data_with_link(als %>%
                     filter(padj < 0.1), "dexseq-sig.csv",
                     "All significant DEXSeq results", dropbox_dir)
```
