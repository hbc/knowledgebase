## Building hybrid murine/transgene reference genome

Heather Wick

This document is based on [prior contributions by James Billingsley](https://github.com/hbc/knowledgebase/blob/master/bcbio/Creating_Hybrid_Mammal_Viral_Reference_Genome.md)

1) <b>Download reference genome from ensembl</b>
  * </b>Example: latest mouse reference (GRCm9) fasta (primary assembly) and gtf downloaded [here](http://useast.ensembl.org/Mus_musculus/Info/Index![image](https://github.com/hbc/knowledgebase/assets/33556230/98d91abd-5cd9-4651-b541-48e3bf413483)
)

2) <b>Acquire transgenes and format to match standard fasta format</b>
  * In this case, genes were provided by the client
  * For our purposes, each transgene was considered to be on its own chromosome, the length of which was the length of the individual gene
  * Editing was done via plain text editor
    
    Format:
      ```
      > [GENE_NAME] dna:chromosom chromosome:[GENOME]:[GENE_NAME]:[CHR_START]:[CHR_END]:1 REF
      BASEPAIRS_ALLCAPS_60_CHARACTERS_WIDE
      ```
    Example:
      ```
      >H2B-GFP dna:chromosome chromosome:GRCm39:H2B-GFP:1:1116:1 REF
      ATGCCAGAGCCAGCGAAGTCTGCTCCCGCCCCGAAAAAGGGCTCCAAGAAGGCGGTGACT
      AAGGCGCAGAAGAAAGGCGGCAAGAAGCGCAAGCGCAGCCGCAAGGAGAGCTATTCCATC
      TATGTGTACAAGGTTCTGAAGCAGGTCCACCCTGACACCGGCATTTCGTCCAAGGCCATG
      GGCATCATGAATTCGTTTGTGAACGACATTTTCGAGCGCATCGCAGGTGAGGCTTCCCGC
      CTGGCGCATTACAACAAGCGCTCGACCATCACCTCCAGGGAGATCCAGACGGCCGTGCGC
      CTGCTGCTGCCTGGGGAGTTGGCCAAGCACGCCGTGTCCGAGGGTACTAAGGCCATCACC
      AAGTACACCAGCGCTAAGGATCCACCGGTCGCCACCATGGTGAGCAAGGGCGAGGAGCTG
      TTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTC
      AGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATC
      TGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGC
      GTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCC
      ATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAG
      ACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGC
      ATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGC
      CACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATC
      CGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCC
      ATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTG
      AGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCC
      GGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA
      ```
      
3) <b>Create GTF for transgenes</b>
  * [Information about GTF file format can be found here](https://useast.ensembl.org/info/website/upload/gff.html)
  * For our purposes, each transgene was considered to be on its own chromosome, the length of which was the length of the individual gene
  * There are only 9 columns. The Attributes in the last column are separated by semicolons/spaces, not tabs.
    
    Format:
      ```
      GENE_NAME	SOURCE	FEATURE	CHR_START	CHR_END	SCORE	STRAND	FRAME	ATTRIBUTES;SEPARATED;BY;SEMI-COLONS;NOT;TABS!;
      ```
    Example:
      ```
      H2B-GFP	unknown	exon	1	1116	.	+	.	gene_id "H2B-GFP"; transcript_id "H2B-GFP"; gene_name "H2B-GFP"; gene_biotype "protein_coding";
      ```
      
4) <b>Concatenate reference fasta and gtf with transgene fasta and gtfs</b>

    Format:
    ```
    cat GENOME.dna.primary_assembly.fa TRANSGENE.fa > GENOME.dna.primary_assembly_TRANSGENE.fa
    cat GENOME.gtf TRANSGENE.gtf > GENOME_TRANSGENE.gtf
    ```
    Example (two transgenes were added):
    ```
    cat Mus_musculus.GRCm39.dna.primary_assembly.fa H2B-GFP.fa tTA.fa > Mus_musculus.GRCm39.dna.primary_assembly_GFP_tTA.fa
    cat Mus_musculus.GRCm39.110.gtf H2B-GFP.gtf tTA.gtf > Mus_musculus.GRCm39.110_GFP_tTA.gtf
    ```
  * <b>Check your formatting!! Sometimes extra new lines or tabs are easy to accidentally add</b>
  
5) <b>Create folder for new reference genome</b>
    ```
    sudo -su bcbio /bin/bash
    cd /n/app/bcbio/1.2.9/genomes/Mmusculus/
    mkdir GRCm39
    ```
  * If you wish to move your new fasta and gtf to the new directory, you may need to move the files to your home directory, then use sudo to sign in as bcbio before copying them to their final location. You may also need to log into an interactive session because the files are quite large
    ```
    cd GRCm39
    srun --pty -p interactive --mem 500M -t 0-06:00
    mv /path/to/home/dir/GENOME.dna.primary_assembly_TRANSGENE.fa .
    mv /path/to/home/dir/GENOME_TRANSGENE.gtf .
    ```
    
6) <b>Run bcbio_setup_genome.py</b>
  * Process is too long to run interactively, so use a script:
    
    Format:
    ```
    #!/bin/bash
    #SBATCH -t 5-00:00              # Runtime in D-HH:MM format
    #SBATCH --job-name=genome_setupbcbio       # Job name
    #SBATCH -c 10                       # cores
    #SBATCH -p medium
    #SBATCH --mem-per-cpu=5G           # Memory needed per CPU or --mem
    #SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
    #SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
    #SBATCH --mail-type=ALL             # Type of email notification (BEGIN, END, FAIL, ALL)
    #SBATCH --mail-user=[USER]k@hsph.harvard.edu
    
    #this script submits bcbio_setup_genome.py
    #must sudo into bcbio first:
    sudo -su bcbio /bin/bash
    
    date
    
    bcbio_setup_genome.py -f GENOME.dna.primary_assembly_TRANSGENE.fa -g GENOME_TRANSGENE.gtf -i bwa star seq -n SPECIES -b GENOME --buildversion BUILD
    
    date
    ```
    Example of actual script can be found here: `/n/app/bcbio/1.2.9/genomes/Mmusculus/GRCm39_GFP_tTA/submit_genome_setup.sh`
    
