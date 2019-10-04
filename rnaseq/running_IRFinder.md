---
title: How to run intron retention analysis
description: This code helps to run IRFinder in the cluster.
category: research
subcategory: rnaseq
tags: [hpc, intro_retention]
---

To run any of these commands, need to activate the bioconda IRFinder environment prior to running script.

1. First script creates reference build required for IRFinder

    ```bash
    #SBATCH -t 24:00:00                    # Runtime in minutes
    #SBATCH -n 4
    #SBATCH -p medium                # Partition (queue) to submit to
    #SBATCH --mem=128G        # 128 GB memory needed (memory PER CORE)
    #SBATCH -o %j.out               # Standard out goes to this file
    #SBATCH -e %j.err               # Standard err goes to this file
    #SBATCH --mail-type=END         # Mail when the job ends

    IRFinder -m BuildRefProcess -r reference_data/
    ```

      >**NOTE:** The files in the `reference_data` folder are sym links to the bcbio ref files and need to be named specifically `genome.fa` and `transcripts.gtf`:
      >
      >`genome.fa -> /n/app/bcbio/biodata/genomes/Hsapiens/hg19/seq/hg19.fa`
      >
      >`transcripts.gtf -> /n/app/bcbio/biodata/genomes/Hsapiens/hg19/rnaseq/ref-transcripts.gtf`

2. Second script (.sh) runs IRFinder and STAR on input file

      ```bash
      #!/bin/bash

      module load star/2.5.4a

      IRFinder -r /path/to/irfinder/reference_data \
      -t 4 -d results \
      $1
      ```

3. Third script (.sh) runs a batch job for each input file in directory

      ```bash
      #!/bin/bash

      for fq in /path/to/*fastq
      do

      sbatch -p medium -t 0-48:00 -n 4 --job-name irfinder --mem=128G -o %j.out -e %j.err --wrap="sh /path/to/irfinder/irfinder_input_file.sh $fq"
      sleep 1 # wait 1 second between each job submission

      done
      ```

4. Fourth script takes output (IRFinder-IR-dir.txt) and uses the replicates to determine differential expression using the Audic and Claverie test (# replicates < 4). analysisWithLowReplicates.pl script comes with the IRFinder github repo clone, so I cloned the repo at https://github.com/williamritchie/IRFinder/. Notes on the Audic and Claverie test can be found at: https://github.com/williamritchie/IRFinder/wiki/Small-Amounts-of-Replicates-via-Audic-and-Claverie-Test.

      ```bash
      #!/bin/bash

      #SBATCH -t 24:00:00                    # Runtime in minutes
      #SBATCH -n 4
      #SBATCH -p medium                # Partition (queue) to submit to
      #SBATCH --mem=128G        # 8 GB memory needed (memory PER CORE)
      #SBATCH -o %j.out               # Standard out goes to this file
      #SBATCH -e %j.err               # Standard err goes to this file
      #SBATCH --mail-type=END         # Mail when the job ends

      analysisWithLowReplicates.pl \
        -A A_ctrl/Pooled/IRFinder-IR-dir.txt A_ctrl/AJ_1/IRFinder-IR-dir.txt A_ctrl/AJ_2/IRFinder-IR-dir.txt A_ctrl/AJ_3/IRFinder-IR-dir.txt \
        -B B_nrde2/Pooled/IRFinder-IR-dir.txt B_nrde2/AJ_4/IRFinder-IR-dir.txt B_nrde2/AJ_5/IRFinder-IR-dir.txt B_nrde2/AJ_6/IRFinder-IR-dir.txt \
        > KD_ctrl-v-nrde2.tab
      ```

5. Output `KD_ctrl-v-nrde2.tab` file can be read directly into R for filtering and results exploration.

6. Rmarkdown workflow (included in report): IRFinder_report.md
