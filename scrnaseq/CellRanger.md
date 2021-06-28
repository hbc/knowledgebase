# When running Cell Ranger on O2
## Shared by Victor
> 10x documentation is the one he uses:

https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger

> For the custom genome:

https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_mr

> For GFP:

https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_mr#marker

## Build custom genome 
> High level steps (based on the 10x tutorial):
1. Download the gtf and fasta files for the species of interest;
2. Filter gtf with `cellranger mkgtf` command; 
3. Create the fasta file for the additional gene (for example, GFP);
4. Create the corresponding gtf file for the additional gene;
5. Append fasta file of the additional gene to the end of the fasta file for the genome;
6. Append gtf file of the additional gene to the end of the gtf file for the genome;
7. Make custom genome with `cellranger mkref` command;

