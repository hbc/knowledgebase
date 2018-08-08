# 3' DGE from LSP demultiplexing example
```
bcl2fastq --adapter-stringency 0.9 --barcode-mismatches 0 --fastq-compression-level 4 --min-log-level INFO --minimum-trimmed-read-length 0 --sample-sheet /n/boslfs/INSTRUMENTS/illumina/180604_NB501677_0276_AHVMT2BGX5/SampleSheet.csv --runfolder-dir /n/boslfs/INSTRUMENTS/illumina/180604_NB501677_0276_AHVMT2BGX5 --output-dir /n/boslfs/ANALYSIS/180604_NB501677_0276_AHVMT2BGX5 --processing-threads 8 --no-lane-splitting --mask-short-adapter-reads 0 --use-bases-mask y*,y*,y*,y*
```