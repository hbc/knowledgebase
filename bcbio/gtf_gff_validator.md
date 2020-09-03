# GTF and GFF valitors
Found during Danesh Moazed consult

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

##GTF validator (perl-based)
- Download: https://mblab.wustl.edu/software.html#evalLink (press 'RELEASES")
- Documentation: https://mblab.wustl.edu/media/software/eval-documentation.pdf

### Usage
```{bash, eval=FALSE}
let's say tar is downloaded and extracted at {/home/eval-2.2.8}
that eval folder is note as {eval}
perl -I {eval} {eval}/validate_gtf.pl {gtf_file} {fasta_file_associated_with_the_gtf}
```
