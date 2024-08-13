# Velocyto documentation

Velocyto is a necessary first step to running RNA velocity

Getting this to work on the server can be a bit fidgety. Things to consider:
* There are multiple wrappers for velocyto, and which one you use depends on:
  * whether you're using data from a simple scRNAseq run or from a multiomics run (even if only using the RNAseq aspect)
  * whether you want to write output to a different folder
* Memory allocation is finnicky with samtools; you need both enough memory for each run and enough memory for each thread, and memory for each thread * number of threads must be less than or equal to memory requested

## Initial set up

**Sources:**
* https://github.com/hbc/tutorials/blob/7ff670c5e3b477da09b6c2e832e05bd43e25448f/scRNAseq/scRNAseq_analysis_tutorial/lessons/velocity.md?plain=1
* https://velocyto.org/velocyto.py/install/index.html

This initial setup uses `virtualenv` rather than `conda`; it is possible to do a conda setup, but personally I ran into issues with that. Mary Piper's existing virtualenv tutorial worked beautifully though, and is recreated below with modifications to reflect the need for updated packages (primarily needing a more updated samtools)

**note: if something is commented out you only need to run it once for initial environment setup. For the actual analysis though, make sure you have the virtual environment loaded with the necessary modules as well**

```
## On command line
# log into interactive node:
    interactive

## Load modules
module load gcc/9.2.0 python/3.8.12 samtools
    
## Create virtual environment
# virtualenv velocyto --system-site-packages
    
# Activate virtual environment
source velocyto/bin/activate
  
# Install tools
# pip3 install numpy scipy cython numba matplotlib scikit-learn h5py click
# pip3 install velocyto
# pip3 install scvelo
```

## Make sure you have proper `10x` gtf and download repeat masking gtf from ensembl

* You can acquire the 10x gtf from the bauer (or other) sequencing center
* the repeat masking gtf can be downloaded from ucsc
  * for example grch38 can be downloaded from https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=611454127_NtvlaW6xBSIRYJEBI0iRDEWisITa&clade=mammal&org=Human&db=0&hgta_group=allTracks&hgta_track=rmsk&hgta_table=rmsk&hgta_regionType=genome&position=&hgta_outputType=gff&hgta_outFileName=hg38_rmsk.gtf

## Setting up your velocyto script -- what wrapper and what memory?

`velocyto run10x` is a very thin wrapper around `veloctyo run` with some important differences and usages:

### Important usage notes for both `velocyto run10x` and `veloctyo run`:
* If you allow it to run samtools to sort the bam files (default) you need write access to the 10x `outs` folder
* Memory requirements for samtools are very particular:
  * In my experience, you MUST provide `--samtools threads` and `--samtools-memory` arguments
	* The `--samtools-memory` is *per node* and the number you specify there times the number of threads you specify in  `--samtools-thread` *must not exceed* the memory requested for the job in `SBATCH --mem`
* You must have write access to 10x `outs` folder in order for samtools to do its sorting, as it needs to be able to write temp files; there is no way to specify a different folder for samtools output, even if you are using a new outfolder for final `.loom` files in `velocyto run`
*   * If need be, copy your 10x `outs` folders to a new directory to get proper permissions
* You can run samtools in advance to presort the files if you want to get around this, but I believe the same memory settings will apply
  * To run: `samtools sort -t CB -O BAM -o cellsorted_possorted_genome_bam.bam possorted_genome_bam.bam`  (see https://velocyto.org/velocyto.py/tutorial/cli.html)
	* Not sure if you would need to regenerate the .bam.bai index or the 10x barcodes file; could not find an answer to this and abandoned pre-sorting as I was running into same issues as with `velocyto run`
	
### `veloctyo run10x`:
* Only the parent folder to the 10x `outs` folder is supplied
* All output will be written to the 10x `outs` folder
* Will NOT work for multiomics 10x output, even if you're using just scRNAseq output from that experiment
  * `run10x` is very particular about how input bam files are named: only accepts scRNAseq bam files named `possorted_genome_bam.bam` (multiomics scRNAseq bam files are named `gex_possorted_genome_bam.bam`) 
	* This peculiarity is NOT clear due to the fact that only the parent folder to the 10x `outs` folder is provided
* Will put output in 10x `outs` folder; cannot specify new out folder

### `veloctyo run`
* You must specify individual input bam files, not just the parent folder to the `outs` folder
  * This allows you to run multiomics files beginning with `gex` as well as regular scRNAseq files beginning with `possorted`
* You must also specify input barcodes file
  * `outs/filtered_feature_bc_matrix/barcodes.tsv.gz`
* you can supply a new outfolder for final `.loom`

## Example code

I ran velocyto on a multi-omics data set. Different samples required different amounts of memory.
For one sample, 160GB total requested memory wasn't enough to process (as it was for previous samples) but increasing to 200GB worked; note I didn't need to increase samtools' memory for this to work.

Here's what the script looked like (please note these directories are real/live so you can go see how the directory is set up -- please don't overwrite them!):

```
#!/bin/sh
#SBATCH -p short
#SBATCH -J velo_day13
#SBATCH -o %x_%j.o
#SBATCH -e %x_%j.e
#SBATCH -t 0-8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --mail-type=ALL

### USAGE:
# sbatch velocyto sample 
# Recommend to run from desired output folder so logs end up where you want them (or change path to log files above)

### Input:
# sample is the name of the parent folder of the outs folder (assumes all paths leading up to folder are the same for all samples)

### Output:
# .loom file in a folder named after $sample

#load before running
#module load gcc/9.2.0 python/3.8.12 samtools
#source ~/velocyto/bin/activate

sample=$1
#copied to local folders
# memory is memory per thread; make sure it's less than memory requested when multiplied by the number of threads!

#everything else
velocyto run -vv -o /n/data1/cores/bcbio/PIs/werner_neuhausser/neuhausser_scRNA-seq_human_embryo_hbc04528/multiome_combined-2023-10-16/RNAvelocity/velocyto_out/${sample}/  \
        -m /n/data1/cores/bcbio/PIs/werner_neuhausser/neuhausser_scRNA-seq_human_embryo_hbc04528/multiome_combined-2023-10-16/RNAvelocity/ref/hg38_rmsk.gtf \
        -b /n/data1/cores/bcbio/PIs/werner_neuhausser/neuhausser_scRNA-seq_human_embryo_hbc04528/multiome_combined-2023-10-16/RNAvelocity/data_copy/${sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
        --samtools-threads 16 \
        --samtools-memory 10000 \
        /n/data1/cores/bcbio/PIs/werner_neuhausser/neuhausser_scRNA-seq_human_embryo_hbc04528/multiome_combined-2023-10-16/RNAvelocity/data_copy/${sample}/outs/gex_possorted_bam.bam \
        /n/data1/cores/bcbio/PIs/werner_neuhausser/neuhausser_scRNA-seq_human_embryo_hbc04528/multiome_combined-2023-10-16/RNAvelocity/ref/genes.gtf
```
