*Single cell analyses require a lot of memory and often fail on the laptops. 
Having R + Seurat installed in a conda environment + interactive session or batch jobs with 50-100G RAM helps.*

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

# 5. Install umap-learn for UMAP clustering
```
pip install umap-learn
```

# 6. Deactivate conda
```
conda deactivate
```

# 7. (Troubleshooting)
- It may ask you to install github token - too many packages loaded from github.
I generated token on my laptop and placed it in ~/.Renviron
- BiocManager::install("slingshot") - I failed to install it due to gsl issues.
- when running a batch job, use source activate r/ source deactivate 
- if conda is trying to write in bcbio cache, check and set cache priority, your home cache should be first:  
`conda info`,  
~/.condarc
```
pkgs_dirs:
    - /home/[UID]/.conda/pkgs
    - /n/app/bcbio/dev/anaconda/pkgs
```
