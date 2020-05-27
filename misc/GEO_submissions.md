Main guide is here:
https://www.ncbi.nlm.nih.gov/geo/info/submission.html
Highh throughput sequencing is here:
https://www.ncbi.nlm.nih.gov/geo/info/seq.html

## Analyst responsibilities
You will need
1) [GEO metadata sheet](https://www.ncbi.nlm.nih.gov/geo/info/seq.html)
2) Raw fastq files
3) Derived files for data  
  a) RNAseq
- raw counts table (as tsv/csv), can put as supplementary file
- TPM (as tsv/csv), can put as supplementary file
- bams are OK too but I have never been asked for them

4) details on the analysis for the GEO metadata sheet, including
- which sequencer was used
- paired or single end reads?
- insert size if paired
- programs and versions used in the analysis, including the bcbio and R portions

Example metadata sheets can be found in this Dropbox folder:
https://www.dropbox.com/sh/88035zd8h9qhvzh/AACmHB7xsXhdgrSyZY42uwLYa?dl=0

For all of the raw and derived data files, you will need to run md5 checksums.

## Researcher responsibilities
The client wil need to give you the details about things that were involved in the experiment and library preparation.
These include
- growth protocol
- treatment protocol 
- extract protocol
- library construction protocol
They will also need to supply the general info about the experiment including:
- title
- summary
- overall design
- who they want to be a contributor

I usually fill out what I can and then send them the metadata sheet with their areas to fill out highlighted.



