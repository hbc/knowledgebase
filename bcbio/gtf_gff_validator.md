# GTF and GFF validators

> ### Use case example: bosTau9 genome 
Building a new genome in bcbio. Reference files were retrieved from NCBI (RefSeq genome and gtf files). There are some additional whitespaces in the file causing errors in the build. Solution: download the gff file instead and validate using `gff3validator` and use as input to bcbio with the added parameter `-gff3`. More info on genometools installs and commands found below. Another option woud be to use the GTF validator (perl-based), also listed below.

> **NOTE:** Initially found during Moazed consult


## GFF Validator
- includes gtf to gff converter
- Download latest from: http://genometools.org/pub/
- Documentation: http://genometools.org/tools.html

### Usage
```{bash, eval=FALSE}
{installed_path}/gt -help
{installed_path}/gt gff3validator {gff_file}
{installed_path}/gt gtf_to_gff3 {gtf_file}
{installed_path}/gt gff3_to_gtf {gff_file}
```
### Errors
Report bugs to https://github.com/genometools/genometools/issues.

## GTF validator (perl-based)
- Download: https://mblab.wustl.edu/software.html#evalLink (press 'RELEASES")
- Documentation: https://mblab.wustl.edu/media/software/eval-documentation.pdf

### Usage
let's say tar is downloaded and extracted at **/home/eval-2.2.8**

That folder is noted as **{eval}** in the code:

```{bash, eval=FALSE}
perl -I {eval} {eval}/validate_gtf.pl -f {gtf_file} {fasta_file_associated_with_the_gtf}
```
'-f' is an option, it creates a fixed file with same title as the origial gtf with '.fixed.gtf' extension.
A custom hg38 gtf ran for an hour.
Memory ran out with 8GB for some reason, so I ran with 64GB just in case. Since it might be using information from the genome.fa extensively.
