# How to run Demuxafy on O2

For detailed instructions and updates on `demuxafy`, see the comprehensive [Read the Docs](https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/index.html#)


## Installation

I originally downloaded the `Demuxafy.sif` singularity image for use on O2 as instructed [here](https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/Installation.html). However, **this singularity image did not pass O2's security checks**. The folks at HMS-RC were kind enough to amend the image for me so that it would pass the security checks. The working image is found at: `/n/app/singularity/containers/aj186/demuxafy.sif` (Victor also has a copy). If you'd like to use it yourself, get in touch with HMS-RC to request it copied on your user folder.

Of note, this singularity image includes a bunch of software, including popscle, demuxlet, freemuxlet, souporcell and other demultiplexing as well as doublet detection tools, so very useful to have installed!


## Input data 

Each tool included in `demuxafy` requires slightly different input (see [Read the Docs](https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/index.html#)). 

For the demultiplexing tools, in most cases, you will need:

- A common SNP genotypes VCF file (pre-processed VCF files can be downloaded [here](https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/DataPrep.html), which is what I did after repeatedly failing to re-generate my own VCF file from the 1000 genome dataset following the provided instructions...)
- A Barcode file (`outs/raw_feature_bc_matrix/barcodes.tsv.gz` from a typical `cellranger count` run)
- A BAM file of aligned single-cell reads (`outs/possorted_genome_bam.bam` from a typical `cellranger count` run)
- Knowledge of the number of samples in the pool you're trying to demultiplex
- Potentially, a FASTA file of the genome your sample was aligned to


## Pre-processing steps

Once you've collated those files, you need to make sure your VCF and BAM files are sorted in the same way. This can be achieved by running the following command after sourcing `sort_vcf_same_as_bam.sh` from the Aerts' Lab popscle helper tool GitHub repo (available [here](https://github.com/aertslab/popscle_helper_tools/blob/master/sort_vcf_same_as_bam.sh)):

```
# Sort VCF file in same order as BAM file
sort_vcf_same_as_bam.sh $BAM $VCF > demuxafy/data/GRCh38_1000G_MAF0.01_ExonFiltered_ChrEncoding_sorted.vcf
```

If you wish to run `freemuxlet` (and possibly other tools I haven't piloted), you will also need to run `dsc-pileup` (available within the singularity image) ahead of `freemuxlet` itself:

```
# Run popscle pileup ahead of freemuxlet
singularity exec $DEMUXAFY popscle dsc-pileup --sam $BAM --vcf demuxafy/data/GRCh38_1000G_MAF0.01_ExonFiltered_ChrEncoding_sorted.vcf --group-list $BARCODES --out $FREEMUXLET_OUTDIR/pileup
``` 

## Workflow

After that, you should be set to run whichever demultiplexing tool you want! I have sample scripts for souporcell, vireo and freemuxlet in the following [GitHub repo](https://github.com/hbc/neuhausser_scRNA-seq_human_embryo_hbc04528/tree/main/pilot_scRNA-seq/demuxafy/scripts)

You also have the option to generate combined results files to contrast results from different software more easily, as described [here](https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/CombineResults.html), and as implemented in the `combine_results.sbatch` script in the GitHub repo linked above.
