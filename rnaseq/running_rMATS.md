## rMATS for differential splicing analysis

* Event-based analysis of splicing (e.g. skipped exon, retained intron, alternative 5' and 3' splice site)
* rMATS handles replicate RNA-Seq data from both paired and unpaired study design
* statistical model of rMATS calculates the P-value and false discovery rate that the difference in the isoform ratio of a gene between two conditions exceeds a given user-defined threshold

Software: https://rnaseq-mats.sourceforge.net/

Paper: https://www.pnas.org/doi/full/10.1073/pnas.1419161111

GitHub: https://github.com/Xinglab/rmats-turbo

### Installation

Issues with the conda build installation provided on the GitHub page  `./build_rmats --conda`. Had problem with shared libraries (" "loading shared libraries" error ). 

Instead install from bioconda. Reference: https://groups.google.com/g/rmats-user-group/c/S1GFEqB9TE8/m/YV9R27CoCwAJ?pli=1

```bash

# Need speicific python version
conda create -n "rMATS_python3.7" python=3.7

conda activate rMATS_python3.7

conda install -c conda-forge -c bioconda rmats=4.1.0

```

If you are running as a script on O2:

```bash
#! /bin/bash

#SBATCH -t 0-24:00       # Runtime
#SBATCH -p medium            # Partition (queue)
#SBATCH -J rmats             # Job name
#SBATCH -o rmats_frag.out             # Standard out
#SBATCH -e rmats_frag.err             # Standard error
#SBATCH --mem=50G     # Memory needed per core
#SBATCH -c 6


# USAGE: For paired-end BAM files;run rMATS

# Define the project path
path=/n/data1/cores/bcbio/PIs/peter_sicinski/sicinski_inhibition_RNAseq_human_hbc04676

# Change directories
cd ${path}/rMATS

# Activate conda env for rmats
source ~/miniconda3/bin/activate
conda init bash
source ~/.bashrc
conda activate rMATS_python3.7 

```
