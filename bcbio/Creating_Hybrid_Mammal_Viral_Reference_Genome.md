## Building hybrid mammal virus genome; human with HHV-1 (HSV-1)
<br>
<br>

https://bcbio-nextgen.readthedocs.io/en/latest/contents/configuration.html#adding-custom-genomes


1.	**Get the genomes** and annotations from ensembl with wget on the transfer node

Genomic sequence: ftp://ftp.ensembl.org/pub/release-102/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz


Annotation: ftp://ftp.ensembl.org/pub/release-102/gff3/homo_sapiens/Homo_sapiens.GRCh38.102.gff3.gz


https://bioinformatics.stackexchange.com/questions/540/what-ensembl-genome-version-should-i-use-for-alignments-e-g-toplevel-fa-vs-p

Use "primary assembly" not "top level."



Client had provided virus gff (for a specific strain, HSV-1_kos) that he had edited to remove repeat sequences, and a similarly edited reference

KOSnorepeat.fa

kosNorepeatFinal.gff


Viral sequences are not available on ensembl (with the exception of SARS-CoV-2)

You can find them at NCBI here:

https://www.ncbi.nlm.nih.gov/nuccore/JQ673480.1

To the right side of the title, there is a "Send to". Click that. Choose "Complete Record" . Under "Choose Destination", choose "File". Under "Format" there is a menu. Choose "GFF3" in that menu and click Create File to get the annotation


<br>
<br>


2.	**Edit the gff file**, to make the human and virus gtf formats match in the concatenated files


First converted the hsv gff3 file (it is gff version 3 in header) to gtf using cufflinks 

use gffread -F -T -o 

		module spider cufflinks
		module load cufflinks/2.2.1
		gffread -h

		gffread  kosNorepeatFinal.gff -F -T -o kosNorepeatFinal.gtf

<br>


This the gtf head

gi|380776962|gb|JQ673480.1|	GenBank	exon	186	860	.	+	.	transcript_id "HHV1gp003.t01"; protein_id "AFE62828.1"; Dbxref "GI:380776965"; Name "UL1"; codon_start "1"; locus_tag "HHV1gp003"; product "envelope glycoprotein L"; translation "length.224";
gi|380776962|gb|JQ673480.1|	GenBank	exon	733	1737	.	+	.	transcript_id "HHV1gp004.t01"; protein_id "AFE62829.1"; Dbxref "GI:380776966"; Name "UL2"; codon_start "1"; locus_tag "HHV1gp004"; product "uracil-DNA glycosylase"; translation "length.334";





Replace the virus fasta header and gtf entries with JQ673480.1, the NCBI HHV-1 kos strain ID


\>JQ673480.1
GGGCCCCCCCCAAAACACACCCCCCGGGGGTCGCGCGCGGCCCTTTAAAGGCGGGCGGCGGGTATATAAA
CCAACGAAAAGCGCGGGAACGGGGATACGGGGCTTGTGTGGCACGACGTCGTGGTTGTGTTACTGGGCAA
ACACTTGGGGACTGTAGGTTTCTGTGGGTGCCGACCCTAGGCGCTATGGGGATTTTGGGTTGGGTCGGGC
TTATTGCCGTTGGGGTTTTGTGTGTGCGGGGGGGCTTGTCTTCAACCGAATATGTTATTCGGAGTCGGGT
GGCTCGAGAGGTGGGGGATATATTAAAGGTGCCTTGTGTGCCGCTCCCGTCTGACGATCTTGATTGGCGT
TACGAGACCCCCTCGGCTATAAACTATGCTTTGATAGACGGTATATTTTTGCGTTATCACTGTCCCGGAT



JQ673480.1	GenBank	exon	186	860	.	+	0	gene_id "HHV1gp003.t01"; transcript_id "HHV1gp003.t01"; protein_id "AFE62828.1"; Dbxref "GI:380776965"; gene_name "UL1"; codon_start "0"; locus_tag "HHV1gp003"; product "envelope glycoprotein L"; translation "length.224"; gbkey "exon,Gene";
JQ673480.1	GenBank	exon	733	1737	.	+	0	gene_id "HHV1gp004.t01"; transcript_id "HHV1gp004.t01"; protein_id "AFE62829.1"; Dbxref "GI:380776966"; gene_name "UL2"; codon_start "0"; locus_tag "HHV1gp004"; product "uracil-DNA glycosylase"; translation "length.334"; gbkey "exon,Gene";



The gtf requires the transcript_id field, and be careful of spacing after the ; in the gtf file

Star will need to see "exon" in the third field to count viral transcripts

<br>
<br>

3.	**gunzip** and **cat** the mammal and virus files 


		cat ensembl_genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa HSV-1_genome/KOSnorepeat.fa > catted_genomes/GRCh38.p13.102_HSV-1_kos.fa


		cat ensembl_genomes/Homo_sapiens.GRCh38.102.gtf HSV-1_genome/kosNorepeatEdited.gtf > catted_genomes/GRCh38.p13.102_HSV-1_kos.gtf

<br>
<br>

4.	**bcbio_setup_genome** to build the genome, it will automatically put it in the right place, here, for example:

/n/shared_db/bcbio/biodata/genomes/Hsapiens/


````
#!/bin/bash

#SBATCH --job-name=setupGenome         # Job name
#SBATCH --partition=priority           # Partition name
#SBATCH --time=1-00:00                 # Runtime in D-HH:MM format
#SBATCH --nodes=1                      # Number of nodes (keep at 1)
#SBATCH --ntasks=1                     # Number of tasks per node (keep at 1)
#SBATCH --cpus-per-task=4              # CPU cores requested per task (change for threaded jobs)
#SBATCH --mem=128G                     # Memory needed per node (total)
#SBATCH --error=jobid_%j.err           # File to which STDERR will be written, including job ID
#SBATCH --output=jobid_%j.out          # File to which STDOUT will be written, including job ID
#SBATCH --mail-type=ALL                # Type of email notification (BEGIN, END, FAIL, ALL)

bcbio_setup_genome.py \
-c 4 \
-f /n/data1/cores/bcbio/PIs/don_coen/Emilia_Vanni/vanni_bulk_RNA-Seq_HSV1_human-mouse_hbc04057/catted_genomes/GRCh38.p13.102_HSV-1_kos.fa \
-g /n/data1/cores/bcbio/PIs/don_coen/Emilia_Vanni/vanni_bulk_RNA-Seq_HSV1_human-mouse_hbc04057/catted_genomes/GRCh38.p13.102_HSV-1_kos.gtf \
-i star seq \
-n Hsapiens \
-b GRCh38.p13.102_HSV-1_kos \
--buildversion GRCh38.p13.102_HSV-1_kos
````
<br>
<br>


5. After running alignment, **view output**, i.e


featureCounts/annotated_combined.counts 

<br>
<br>

6.	**filter** out multimapped reads, **sort**, re-run featureCounts

use ready.bam files



		bamtools filter -tag NH:1 -in my.bam -out my.filtered.bam

		bamtools filter -tag NH:1 -in LIB042801-ready.bam -out post_filt_LIB042801_ready.bam


		samtools sort post_filt_LIB042801_ready.bam -o sample.sorted.LIB042801_ready.bam


<br>


or do them all

cp all the ready.bam to a dir

filter and sort


````
#!/bin/bash

#SBATCH --job-name=bamtoolsfilter      # Job name
#SBATCH --partition=priority           # Partition name
#SBATCH --time=1-00:00                 # Runtime in D-HH:MM format
#SBATCH --nodes=1                      # Number of nodes (keep at 1)
#SBATCH --ntasks=1                     # Number of tasks per node (keep at 1)
#SBATCH --cpus-per-task=4              # CPU cores requested per task (change for threaded jobs)
#SBATCH --mem=128G                     # Memory needed per node (total)
#SBATCH --error=jobid_%j.err           # File to which STDERR will be written, including job ID
#SBATCH --output=jobid_%j.out          # File to which STDOUT will be written, including job ID
#SBATCH --mail-type=ALL                # Type of email notification (BEGIN, END, FAIL, ALL)

module load gcc/6.2.0

module load bamtools/2.4.1

module load samtools/1.10


for file in *.bam; do
    bamtools filter -tag NH:1 -in $file -out "`basename $file .bam`filtered.bam"
done


for file in *filtered.bam; do

    samtools sort $file -o "`basename $file .bam`sorted.bam"
done
````
<br>
<br>

7. run **featureCounts**


````
#!/bin/bash

#SBATCH --job-name=featureCounts       # Job name
#SBATCH --partition=priority           # Partition name
#SBATCH --time=1-00:00                 # Runtime in D-HH:MM format
#SBATCH --nodes=1                      # Number of nodes (keep at 1)
#SBATCH --ntasks=1                     # Number of tasks per node (keep at 1)
#SBATCH --cpus-per-task=4              # CPU cores requested per task (change for threaded jobs)
#SBATCH --mem=128G                     # Memory needed per node (total)
#SBATCH --error=jobid_%j.err           # File to which STDERR will be written, including job ID
#SBATCH --output=jobid_%j.out          # File to which STDOUT will be written, including job ID
#SBATCH --mail-type=ALL                # Type of email notification (BEGIN, END, FAIL, ALL)

featureCounts -T 4 -p -t exon -g gene_id -F GTF \
  -a /n/data1/cores/bcbio/PIs/don_coen/Emilia_Vanni/vanni_bulk_RNA-Seq_HSV1_human-mouse_hbc04057/catted_genomes/GRCh38.p13.102_HSV-1_kos.gtf \
  -o featurecounts \
  LIB042801-readyfilteredsorted.bam LIB042802-readyfilteredsorted.bam LIB042803-readyfilteredsorted.bam
````





then **cut** to get output

		cut -f1,7-8 featurecounts.txt > cutFeaturecounts.txt



