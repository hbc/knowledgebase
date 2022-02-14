## Hybrid Assembly Strategies (smaller/prokaryotic genomes)


## Hybrid Assembly Strategies (larger/eukaryotic genomes)

### Oxford nanopore and Illumina -

> Note that these were suggestions for an **algal assembly** from Dr. Chris Fields and Kim Walden at [UIUC's bioinformatics core, HPCBio](https://hpcbio.illinois.edu/).

* Get a good estimate of the genome size using your illumina reads and [Genome Scope](http://qb.cshl.edu/genomescope/). You will need to get a kmer histogram from your illumina data to use as input to Genome Scope, and you can use [KMC](https://github.com/refresh-bio/KMC) or [Jellyfish](http://www.genome.umd.edu/jellyfish.html) for that.
* Workflow
  * Assemble Nanopore reads using something like [wtdbg2](https://github.com/ruanjue/wtdbg2) or [Flye](https://github.com/fenderglass/Flye) (might need to use multiple assemblers and test which works best)
  * Do a first round of assembly polishing with [nanopolish](https://github.com/jts/nanopolish), using only nanopore data
  * Do a second round of assembly polishing with [Racon](https://github.com/isovic/racon) or [Pilon](https://github.com/broadinstitute/pilon/wiki), using illumina data this time
* Use [BUSCO](https://busco.ezlab.org/) for assessment
* If your Nanopore data was generated using tech from before 2020 using older flow cells, you may have to re run the basecalling using something like [Guppy](https://esr-nz.github.io/gpu_basecalling_testing/gpu_benchmarking.html) (speedy if you can use GPUs)

## Genome Annotation tools

* MAKER
* Braker
* GeneMark
* antiSMASH
* PGAP (prokaryotes)

## Assembly assessment

* [BUSCO](https://busco.ezlab.org/)
