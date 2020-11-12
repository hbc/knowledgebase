# Bcbio workflow by Mary Piper
https://github.com/marypiper/bcbio_rnaseq_workflow/blob/master/bcbio_rna-seq_workflow.md

**Documentation for bcbio:** [bcbio-nextgen readthedocs](http://bcbio-nextgen.readthedocs.org/en/latest/contents/pipelines.html#rna-seq)

## Set-up
1. Follow instructions for starting an analysis using https://github.com/hbc/knowledgebase/blob/master/admin/setting_up_an_analysis_guidelines.md.

3. Download fastq files from facility to data folder
	
	- Download fastq files from a non-password protected url
		- `wget --mirror url` (for each file of sample in each lane)
   	 	- Rory's code to concatenate files for the same samples on multiple lanes: 
    
    			barcodes="BC1 BC2 BC3 BC4"
    			for barcode in $barcodes
    			do
    			find folder -name $barcode_*R1.fastq.gz -exec cat {} \; > data/${barcode}_R1.fastq.gz
    			find folder -name $barcode_*R2.fastq.gz -exec cat {} \; > data/${barcode}_R2.fastq.gz
    			done

   	- Download from password protected FTP such as Dana Farber
		- `wget -r <FTP address of folder> --user <username>  --password <pwd> <destination>`
	
	- Download fastq files from BioPolymers: 
   		- `rsync -avr username@bpfngs.med.harvard.edu:./folder_name .`
   		
   		--OR--
   		
		- `sftp username@bpfngs.med.harvard.edu`
		- `cd` to correct folder
		- `mget *.tab`
		- `mget *.bz2`
		
	- Download from the Broad using Aspera:
		- To download data I use this [script](https://github.com/marypiper/bcbio_rnaseq_workflow/blob/master/aspera_connect_lsf).

4. Create metadata in Excel create sym links by concatenate("ln -s ", column $A2 with path_to_where_files_are_stored, " ", column with name of sym link $D2). Can extract parts of column using delimiters in Data tab column to text.

5. Save Excel as text and replace ^M with new lines in vim:

	`:%s/<Ctrl-V><Ctrl-M>/\r/g`

6. Settings for bcbio- make sure you have following settings in `~/.bashrc` file:
 
 ```bash
    unset PYTHONHOME
    unset PYTHONPATH
    export PATH=/n/app/bcbio/tools/bin:$PATH
 ```
    
7. Within the `meta` folder, add your comma-separated metadata file (`projectname_rnaseq.csv`)
	- first column is `samplename` and is the names of the fastq files as they appear in the directory (should be the file name without the extension (no .fastq or R#.fastq for paired-end reads))
	- second column is `description` and is unique names to call samples - provide the names you want to have the samples called by 
	- **FOR CHIP-SEQ** need additional columns:
		- `phenotype`: `chip` or `input` for each sample
		- `batch`: batch1, batch2, batch3, ... for grouping each input with it's appropriate chip(s)
	- additional specifics regarding the metadata file: [http://bcbio-nextgen.readthedocs.org/en/latest/contents/configuration.html#automated-sample-configuration](http://bcbio-nextgen.readthedocs.org/en/latest/contents/configuration.html#automated-sample-configuration) 
        
8. Within the `config` folder, add your custom Illumina template
    - Example template for human RNA-seq using Illumina prepared samples (genome_build for mouse = mm10, human = hg19 or hg38 (need to change star to hisat2 if using hg38):

	```yaml
	# Template for mouse RNA-seq using Illumina prepared samples
	---
	details:
	  - analysis: RNA-seq
	    genome_build: mm10
	    algorithm:
	      aligner: star
	      quality_format: standard
	      strandedness: firststrand
	      tools_on: bcbiornaseq
	      bcbiornaseq:
		organism: mus musculus
		interesting_groups: [genotype]
	upload:
	  dir: /n/data1/cores/bcbio/PIs/vamsi_mootha/hbc_mootha_rnaseq_of_metabolite_transporter_KO_mouse_livers_hbc03618_1/bcbio_final
	```

	- List of genomes available can be found by running `bcbio_setup_genome.py`
	- strandedness options: `unstranded`, `firststrand`, `secondstrand`
	- Additional parameters can be found: [http://bcbio-nextgen.readthedocs.org/en/latest/contents/configuration.html#automated-sample-configuration](http://bcbio-nextgen.readthedocs.org/en/latest/contents/configuration.html#automated-sample-configuration) 
	- Best practice templates can be found: [https://github.com/chapmanb/bcbio-nextgen/tree/master/config/templates](https://github.com/chapmanb/bcbio-nextgen/tree/master/config/templates)

 
9. Within the `data` folder, add all your fastq files to analyze.

## Analysis

1. Go to `/n/scratch2/your_ECommonsID/PI` and create an `analysis` folder. Change directories to `analysis` folder and create the full Illumina instructions using the Illumina template created in Set-up: step #6.
    - `srun --pty -p interactive -t 0-12:00 --mem 8G bash` start interactive job
    - `cd path-to-folder/analysis` change directories to analysis folder
    - `bcbio_nextgen.py -w template /n/data1/cores/bcbio/PIs/path_to_templates/star-illumina-rnaseq.yaml /n/data1/cores/bcbio/PIs/path_to_meta/*-rnaseq.csv /n/data1/cores/bcbio/PIs/path_to_data/*fastq.gz` run command to create the full yaml file

2. Create script for running the job (in analysis folder)

For a larger job:

	```bash
	#!/bin/sh
	#SBATCH -p medium
	#SBATCH -J mootha
	#SBATCH -o run.o
	#SBATCH -e run.e
	#SBATCH -t 0-100:00
	#SBATCH --cpus-per-task=1
	#SBATCH --mem-per-cpu=8G
	#SBATCH --mail-type=ALL
	#SBATCH --mail-user=piper@hsph.harvard.edu
	
	export PATH=/n/app/bcbio/tools/bin:$PATH
	
	/n/app/bcbio/dev/anaconda/bin/bcbio_nextgen.py ../config/\*\_rnaseq.yaml -n 48 -t ipython -s slurm -q medium -r t=0-100:00 --timeout 300 --retries 3
	```
	
For a smaller job, it might be faster in overall time to just run the job on the priority queue. If you only have a few samples, and your fairshare score is low, running on the priority queue could end up being faster since you will quickly get a job there and not have to wait.

```bash
#!/bin/sh
#SBATCH -p priority
#SBATCH -J mootha
#SBATCH -o run.o
#SBATCH -e run.e
#SBATCH -t 0-100:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=64G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=piper@hsph.harvard.edu
export PATH=/n/app/bcbio/tools/bin:$PATH
/n/app/bcbio/dev/anaconda/bin/bcbio_nextgen.py ../config/\*\_rnaseq.yaml -n 8
```

3. Go to work folder and start the job - make sure in an interactive session 

	```bash
	cd /n/scratch2/path_to_folder/analysis/\*\_rnaseq/work
	sbatch ../../runJob-\*\_rnaseq.slurm
	```

### Exploration of region of interest

1. The bam files will be located here: `path-to-folder/*-rnaseq/analysis/*-rnaseq/work/align/SAMPLENAME/NAME_*-rnaseq_star/` # needs to be updated

2. Extracting interesting region (example)
	- `samtools view -h -b  sample1.bam "chr2:176927474-177089906" > sample1_hox.bam`

	- `samtools index sample1_hox.bam`


## Mounting bcbio

`sshfs mp298@transfer.orchestra.med.harvard.edu:/n/data1/cores/bcbio ~/bcbio -o volname=bcbio -o follow_symlinks`
