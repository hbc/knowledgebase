# How to run Demuxafy on O2

For detailed instructions and updates on `demuxafy`, see the comprehensive [Read the Docs](https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/index.html#)


## Installation

I originally downloaded the `Demuxafy.sif` singularity image for use on O2 as instructed [here](https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/Installation.html). However, **this singularity image did not pass O2's security checks**. The folks at HMS-RC were kind enough to amend the image for me so that it would pass the security checks. The working image is found at: `/n/app/singularity/containers/Demuxafy.sif` allowing anyone to use it.

Of note, this singularity image includes a bunch of software, including popscle, demuxlet, freemuxlet, souporcell and other demultiplexing as well as doublet detection tools, so very useful to have installed!


## Input data 

Each tool included in `demuxafy` requires slightly different input (see [Read the Docs](https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/index.html#)). 

For the demultiplexing tools, in most cases, you will need:

- A common SNP genotypes VCF file (pre-processed VCF files can be downloaded [here](https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/DataPrep.html), which is what I did after repeatedly failing to re-generate my own VCF file from the 1000 genome dataset following the provided instructions...)
- A Barcode file (`outs/raw_feature_bc_matrix/barcodes.tsv.gz` from a typical `cellranger count` run)
- A BAM file of aligned single-cell reads (`outs/possorted_genome_bam.bam` from a typical `cellranger count` run)
- Knowledge of the number of samples in the pool you're trying to demultiplex
- Potentially, a FASTA file of the genome your sample was aligned to

_NOTE_: When working from a multiplexed dataset (e.g. cell hashing experiment), you may have to re-run `cellranger count` instead of `cellranger multi` to generate the proper barcodes and BAM files. In addition, it may be necessary to use the `barcodes.tsv.gz` file from the `filtered_feature_bc_matrix` (instead of raw) in such cases (see for example this [issue](https://github.com/wheaton5/souporcell/issues/128) when running `souporcell`). 


## Pre-processing steps

Once you've collated those files, you need to make sure your VCF and BAM files are sorted in the same way. This can be achieved by running the following command after sourcing `sort_vcf_same_as_bam.sh` from the Aerts' Lab popscle helper tool GitHub repo (available [here](https://github.com/aertslab/popscle_helper_tools/blob/master/sort_vcf_same_as_bam.sh)):

```
# Sort VCF file in same order as BAM file
sort_vcf_same_as_bam.sh $BAM $VCF > demuxafy/data/GRCh38_1000G_MAF0.01_ExonFiltered_ChrEncoding_sorted.vcf
```

### dsc pileup

If you wish to run `freemuxlet` (and possibly other tools I haven't piloted), you will also need to run `dsc-pileup` (available within the singularity image) ahead of `freemuxlet` itself. For larger samples (>30k cells), it also helps (= significantly speeds up computational time, from several days to a couple of hours) to pre-filter the BAM file using another of the Aerts' Lab popscle helper tool scripts: `filter_bam_file_for_popscle_dsc_pileup.sh` (available [here](https://github.com/aertslab/popscle_helper_tools/blob/master/filter_bam_file_for_popscle_dsc_pileup.sh))

```
# [OPTIONAL but recommended]
module load gcc/9.2.0 samtools/1.14
scripts/filter_bam_file_for_popscle_dsc_pileup.sh $BAM $BARCODES demuxafy/data/GRCh38_1000G_MAF0.01_ExonFiltered_ChrEncoding_sorted.vcf demuxafy/data/possorted_genome_bam_filtered.bam

# Run popscle pileup ahead of freemuxlet
singularity exec $DEMUXAFY popscle dsc-pileup --sam demuxafy/data/possorted_genome_bam_filtered.bam --vcf demuxafy/data/GRCh38_1000G_MAF0.01_ExonFiltered_ChrEncoding_sorted.vcf --group-list $BARCODES --out $FREEMUXLET_OUTDIR/pileup
``` 

_NOTE_: When running the dsc-pileup step on O2, at some point the job might get stalled despite no error message being issued. From my experience, this usually means that the requested memory needs to be increased (I used 48G-56G for most samples I processed, and encountered issues when lowering down to 32G). After filtering the BAM file and with the appropriate amount of memory available, the dsc-pileup step usually completes within 2-3 hours.


## Workflow

After that, you should be set to run whichever demultiplexing tool you want! See sample scripts for a simple case (small 10X study) in the following [GitHub repo](https://github.com/hbc/neuhausser_scRNA-seq_human_embryo_hbc04528/tree/main/pilot_scRNA-seq/demuxafy/scripts); and for a more complex case (large study, multiplexed 10X data using cell hashing) [here](https://github.com/hbc/hbc_10xCITESeq_Pregizer-Visterra-_hbc04485/tree/main/demuxafy/scripts)

You also have the option to generate combined results files to contrast results from different software more easily, as described [here](https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/CombineResults.html), and as implemented in the `combine_results.sbatch` script in the first GitHub repo linked above.
