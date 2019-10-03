# Counting cells with bcbio for inDrops3 data - short version

## 1. Check reference genome and transcriptome
(O2 or another bcbio installation):
- mm10 reference genome: /n/shared_db/bcbio/biodata/genomes/Mmusculus/mm10
- transcriptome_fasta: /n/shared_db/bcbio/biodata/genomes/Mmusculus/mm10/rnaseq/ref-transcripts.fa
- transcriptome_gtf: /n/shared_db/bcbio/biodata/genomes/Mmusculus/mm10/rnaseq/ref-transcripts.gtf

## 2. Create bcbio project structure in your /scratch
```
mkdir sc_mouse
cd sc_mouse
mkdir config input final work
```

## 3. Copy(move) input data (fastq) to sc_mouse/input
```
KM_1.fq.gz
KM_2.fq.gz
KM_3.fq.gz
KM_4.fq.gz
```

## 4. Create `sc_mouse/config/sample_barcodes.csv`
```
TCTCTCCG,S01
GCGTAAGA,S02
CCTAGAGT,S03
TCGACTAG,S04
TTCTAGAG,S05
```

## 5. Create `sc_mouse/config/sc-mouse.yaml`
```
details:
- algorithm:
    cellular_barcode_correction: 1
    minimum_barcode_depth: 1000
    sample_barcodes: /full/path/sc_mouse/config/sample_barcodes.csv
    transcriptome_fasta: /n/shared_db/bcbio/biodata/genomes/Mmusculus/mm10/rnaseq/ref-transcripts.fa
    transcriptome_gtf: /n/shared_db/bcbio/biodata/genomes/Mmusculus/mm10/rnaseq/ref-transcripts.gtf
    umi_type: harvard-indrop-v3
  analysis: scRNA-seq
  description: PI_name
  files:
  - /full/path/sc_mouse/input/KM_1.fq.gz
  - /full/path/sc_mouse/input/KM_2.fq.gz
  - /full/path/sc_mouse/input/KM_3.fq.gz
  - /full/path/sc_mouse/input/KM_4.fq.gz
  genome_build: mm10
  metadata: {}
fc_name: sc-mouse
upload:
  dir: /full/path/sc_mouse/final
```

## 6. Create `sc_mouse/config/bcbio.sh`
```
#!/bin/bash

# https://slurm.schedmd.com/sbatch.html

#SBATCH --partition=priority        # Partition (queue)
#SBATCH --time=5-00:00              # Runtime in D-HH:MM format
#SBATCH --job-name=km            # Job name
#SBATCH -c 20
#SBATCH --mem-per-cpu=5G            # Memory needed per CPU
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL             # Type of email notification (BEGIN, END, FAIL, ALL)

bcbio_nextgen.py ../config/sc-mouse.yaml -n 20
```

## 7. Run bcbio
```
cd sc_mouse_work
sbatch ../config/bcbio.sh
```

## 1a (Optional). 
If you care, download fresh transcriptome annotation from Gencode (https://www.gencodegenes.org/mouse/)
(it has chrom names with chr matching mm10 assembly).
```
cd sc_mouse/input
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.annotation.gff3.gz
gunzip gencode.vM23.annotation.gtf.gz
gffread -g /n/shared_db/bcbio/biodata/genomes/Mmusculus/mm10/seq/mm10.fa gencode.vM23.annotation.gtf -x gencode.vM23.annotation.cds.fa
```
update sc_mouse/config/sc_mouse.yaml:
```
transcriptome_fasta: gencode.vM23.annotation.cds.fa
transcriptome_gtf: gencode.vM23.annotation.gtf
```

## 3a (Optional). 
Merge multiple lanes for every read:
```
zcat KM_lane1_R1.fastq KM_lane2_R1.fastq.gz | gzip > KM_1.fq.gz
...
```

## References
- [Even shorter guide](https://github.com/bcbio/bcbio-nextgen/blob/master/config/templates/indrop-singlecell.yaml)
- [Much more comprehensive guide](https://github.com/hbc/tutorials/blob/master/scRNAseq/scRNAseq_analysis_tutorial/lessons/01_bcbio_run.md)
