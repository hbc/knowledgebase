# 3' DGE from LSP demultiplexing example
```
bcl2fastq --adapter-stringency 0.9 --barcode-mismatches 0 --fastq-compression-level 4 --min-log-level INFO --minimum-trimmed-read-length 0 --sample-sheet /n/boslfs/INSTRUMENTS/illumina/180604_NB501677_0276_AHVMT2BGX5/SampleSheet.csv --runfolder-dir /n/boslfs/INSTRUMENTS/illumina/180604_NB501677_0276_AHVMT2BGX5 --output-dir /n/boslfs/ANALYSIS/180604_NB501677_0276_AHVMT2BGX5 --processing-threads 8 --no-lane-splitting --mask-short-adapter-reads 0 --use-bases-mask y*,y*,y*,y*
```

# Illumina instrument by FASTQ read name
- @HWI-Mxxxx or @Mxxxx - MiSeq
- @Kxxxx - HiSeq 3000(?)/4000
- @Nxxxx - NextSeq 500/550
- @Axxxxx - NovaSeq
- @HWI-Dxxxx - HiSeq 2000/2500
- AAXX, @HWUSI - GAIIx
- BCXX = HiSeq v1.5
- ACXX = HiSeq High-Output v3
- ANXX = HiSeq High-Output v4
- ADXX = HiSeq RR v1
- AMXX, BCXX =HiSeq RR v2
- ALXX = HiSeqX
- BGXX, AGXX = High-Output NextSeq
- AFXX = Mid-Output NextSeq

# Illumina BaseSpace CLI
https://developer.basespace.illumina.com/docs/content/documentation/cli/cli-overview
