---
title: Single Cell Installation with conda on O2
description: This code helps with installing tools for Single Cell analysis with conda.
category: research
subcategory: scrnaseq
tags: [hpc, R, conda]
---

# 1. Use conda from bcbio
```
which conda
/n/app/bcbio/dev/anaconda/bin/conda
conda --version
conda 4.6.14
```

# 2. Create and setup r conda environment
```
conda create -n r r-essentials r-base zlib pandoc
conda init bash
conda config --set auto_activate_base false
. ~/.bashrc 
```

# 3. Activate conda env
```
conda activate r
which R

```

# 4. Install packages from within R

## 4.1 Install Seurat
```
R
install.packages("Seurat")
library(Seurat)
q()
```

## 4.2 Install Monocle
```
R
install.packages(c("BiocManager", "remotes"))
BiocManager::install("monocle")
q()
```

## 4.3 Install liger
```
R
install.packages("devtools")
library(devtools)
install_github("MacoskoLab/liger")
library(liger)
q()
```

# 5. Deactivate conda
```
conda deactivate
```

# 6. (Troubleshooting)
- It may ask you to install github token - too many packages loaded from github.
I generated token on my laptop and placed it in ~/.Renviron
- BiocManager::install("slingshot") - I failed to install it due to gsl issues.
