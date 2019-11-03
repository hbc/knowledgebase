# Setting up an analysis Guidelines
## Initial folder and git setup on on O2/Odyssey
1. Setup repo on code.harvard.edu (unless PI wants it to be public, then setup repo on Github). Use same name explained before with hbc_ prefix:
`hbc_$technology_of_$pilastname_$intervention_on_$tissue_in_$organism_$hbccode`
2. If not present already, make a folder on the server in the PIs directory using this format: `$pifirstname_pilastname`
3. Clone the repo inside this directory
4. Go inside the repo directory, and setup subfolders called:  

	#### data (for raw data)

	#### meta (for extra, unformatted, sample metadata)

	#### templates (bcbio config files)
	#### docs (other information that you might want to keep near the data)


	## Notes on folders
	- these are the FASTQ or similar file types for other technologies.
	- In the data folder, you can have the actual downloaded files, or symlinks to someone else’s downloaded files

	### meta
	- munge the metadata into the format for bcbio and give it a name that will tell you about the particular bcbio run.
	- Here is a simple example of a metadata file with only replicate and genotype as covariates:
https://docs.google.com/spreadsheets/d/18h6qPc7_rGzyg2gTbgyg5Nmo00zBikXBJGzDl9QRmXY/edit?usp=sharing
	- The stem of this file (filename without extension) will be used to name the folder with your final bcbio results
	- Typically, I just call it “bcbio.csv”
	- if I was to do another bcbio run with a different genome than before (Flybase for example), I give it a new, descriptive name (eg. “bcbio_flybase.csv”)


	### templates
	- These are the bcbio templates.
	- you can grab one from a previous project or download from the bcbio repo
	- you can also just run `bcbio_download_template rnaseq` (for example) to get the template for your particular technology
	- The final term in the command will be used to match against the following templates and any having any overlap will be downloaded:		
		~~~
		freebayes-variant.yaml   
		illumina-chipseq.yaml   
		illumina-rnaseq.yaml     
		indrop-singlecell.yaml
		tumor-paired.yaml
		gatk-variant.yaml       
		illumina-fastrnaseq.yaml
		illumina-srnaseq.yaml    
		noalign-variant.yaml
		~~~

		So for example, `bcbio_download_template gatk` will only pull down the gatk-variant.yaml template, but `bcbio_download_template ill` will pull down every template that has illlumina in its’ name.

	- Modify the template file as needed to fit the needs of your analysis. If running bcbio_o2 below, there is no need to modify the `upload:dir` variable

5. Setup .gitignore file to ignore all files you don’t want to sync. In theory we only use git to store code and very small files. (Ignore bcb final folder, and data folder) [*FUTURE: NEED EXAMPLES*]*

## Running bcbio


At the end of the run, you will have a directory structure that looks something like this:

~~~   
├── Homo_sapiens.GRCh38.92.gtf     
├── Homo_sapiens.GRCh38.92-tx2gene.tsv     
├── Homo_sapiens.GRCh38.cdna.all.fa    
├── indrop-rnaseq.yaml   
├── metadata   
│   ├── lane1_NoIndex_L001   
│   ├── lane2_NoIndex_L002      
├── sc-human   
│   ├── config   
│   │   ├── sscc-human.csv   
│   │   ├── sc-human-template.yaml   
│   │   ├── sc-human.yaml   
│   │   └── sc-human.yaml.bak2018-07-12-14-38-13   
│   └── final   
│       ├── 2018-07-19_sc-human   
│       │   ├── bcbio-nextgen-commands.log   
│       │   ├── bcbio-nextgen.log   
│       │   ├── bcb.rds   
│       │   ├── data_versions.csv   
│       │   ├── metadata.csv   
│       │   ├── programs.txt   
│       │   ├── project-summary.yaml   
│       │   ├── tagcounts-dupes.mtx   
│       │   ├── tagcounts-dupes.mtx.colnames   
│       │   ├── tagcounts-dupes.mtx.rownames   
│       │   ├── tagcounts.mtx   
│       │   ├── tagcounts.mtx.colnames   
│       │   └── tagcounts.mtx.rownames   
│       ├── lane1-AGCTTTCT   
│       │   ├── lane1-AGCTTTCT-barcodes-filtered.tsv   
│       │   ├── lane1-AGCTTTCT-barcodes.tsv   
│       │   ├── lane1-AGCTTTCT.mtx   
│       │   ├── lane1-AGCTTTCT.mtx.colnames   
│       │   ├── lane1-AGCTTTCT.mtx.rownames   
│       │   └── lane1-AGCTTTCT-transcriptome.bam   
│       ├── lane2-AAGAGCGT   
│       │   ├── lane2-AAGAGCGT-barcodes-filtered.tsv   
│       │   ├── lane2-AAGAGCGT-barcodes.tsv   
│       │   ├── lane2-AAGAGCGT.mtx   
│       │   ├── lane2-AAGAGCGT.mtx.colnames   
│       │   ├── lane2-AAGAGCGT.mtx.rownames   
│       │   └── lane2-AAGAGCGT-transcriptome.bam      
│       └── mtx.tar.gz   
└── sc-human.csv  
~~~


- Copy project template yaml file to the main directory of consult on server
- Copy a slurm batch file for your analysis to the main directory of consult on server
- Edit both to reflect the properties of your consult, making sure to have the final upload dir point to an appropriate folder

- Here’s an example of an RNAseq template (project-template.yaml) for bcbio:

#### Template for human RNA-seq using Illumina prepared samples
~~~
details:      
  - analysis: RNA-seq      
		genome_build: BDGP6      
	    algorithm:      
     	aligner: star      
	    quality_format: Standard      
       	trim_reads: False      
       	adapters: [truseq, polya]      
       	strandedness: unstranded
	upload:
		dir: /n/data1/cores/bcbio/PIs/mel_feany/RNAseq_of_different_genotypes_in_Drosophila_brain/bcbio/final
  ~~~

- And here is an example slurm batch file (run_slurm.sh):
```
#!/bin/sh
#SBATCH -p short
#SBATCH -J feany
#SBATCH -o run.o
#SBATCH -e run.e
#SBATCH -t 0-12:00
#SBATC H --cpus-per-task=1
#SBATCH --mem=8000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jnhutchinson@gmail.com
export PATH=/n/app/bcbio/tools/bin:$PATH
/n/app/bcbio/dev/anaconda/bin/bcbio_nextgen.py ../config/bcbio_ensembl.yaml -n 72 -t ipython -s slurm -q short -r t=0-6000:00 --tag feany
```

####  Setting up on scratch

It's a good idea to get bcbio to use a work directory that is on scratch. One way to do this is to run bcbio's templating script and then replace the work folder with a symlink to a folder on scratch. keep your work direcotyr  variant has the nice feature of automatically putting the working directory on the scratch drive, so that our storage space doesn’t blow up

Another more convoluted route is just to run everything on scratch and copy over the final files:

- Run the actual bcbio run on Scratch and only upload the final bcbio results to the main directories
- Setup directories PI folder and consult subfolder on Scratch the same as you did in steps 1 and 2 of the Raw Data steps
- Symlink to the raw data directory
- Symlink to the bcbio metadata file you made in step 2 of Metadata section
- Symlink to the config files you made
- Generate a the template for the run on scrath using bcbio’s template options
eg. bcbio_nextgen.py -w template project-template.yaml bcbiol.csv ./data/*gz
- Go into bcbio/config folder, check to make sure the bcbio.yaml file has the proper info (i.e. sample names are correct, all samples are there)
- Go into the bcbio/work folder and symlink the slurm batch file
- Run the slurm batch file eg. sbatch run_slurm.sh
- Monitor job progress until complete


## Running an Analysis

### Create the bcbioRNASeq object
- Decide where you are going to create the bcbio object
	- If possible, create the bcbioRNASeq object on the server so you have a backup of it. You can run bcbioRNASeq on the cluster using server modules or any other methods, but the easiest is to do it via conda  (for conda installation, see bcbioRNASeq README.md).
	- It can also be run in your local computer if you mount O2 filesystem using sshfs. If you do this, it’s a good idea to sync the object up to the server after step 2.
- Load bcbio results into bcbioRNASeq object and save out to the dated bcbio final project folder (`date_projectname`) (this folder will be the dated subfolder in the final subfolder which is  a subfolder of a folder named after the bcbio metadata file you made above, for a better explanation, try looking at the directory structure above)
- bcbio can now do this automatically with the proper setup in the yaml file...for instance the follwoing tempalte will run generate a bcbioRNASeq object and run the QC template on it with the interesting groups of day and genotype highlighted:

```
   # Template for mouse RNA-seq using Illumina prepared samples
   ---
   details:
     - analysis: RNA-seq
       genome_build: mm10
       algorithm:
         aligner: star
         quality_format: standard
         strandedness: unstranded
        tools_on: bcbiornaseq
        bcbiornaseq:
          organism: mus musculus
          interesting_groups: [day, genotype]
  upload:
    dir: ../bcbio_final
```


###   Sync git repo to local computer
- Clone the repo to your local drive (not in a Dropbox or Google Drive folder!). For consitency, you may want to replicate the server folder structure by putting it inside a folder named after the PI (i.e. `$pifirsthane_$pilastname`)
- Create a README with any information needed to repeat the analysis. For instance, complex consults where you have to make custom scripts/database preparation etc.
- Commit README and push to repo.
- Make a “report” subfolder inside it
- Setup an Rstudio project inside the report folder
- Load up a QC template for bcbioRNASeq and save in report folder as something like QC.Rmd
- Commit QC code to git, and push to repo
- Do analysis and update git as you go      



## Sharing results

- we mainly share our results with researchers on Dropbox
- the Dropbox folders have the same naming scheme as the server folders within the PI directory. Create a folder with the same name as the subfolder under the PIs folder.

You can share via Dropbox by either:

1)  using code within `r2dropSmart::sync` function (install from github lpantano/r2dropSmart)

- do the next two steps once per project directory

	~~~
	library(rdrop2)        
	drop_auth()
	~~~
- this will launch your browser and request access to your Dropbox account. You will be prompted to log in if you aren't already logged in
- once completed, close your browser window and return to R to complete authentication
the credentials are automatically cached (you can prevent this if you’d like, see the rdrop2 do) for future use

- If you wish to save the tokens, for local/remote use:

	~~~
	token <- drop_auth()
	saveRDS(token, file = "token.rds")
	~~~

- To sync your results

	~~~
	library(r2dropSmart)
	token <- readRDS("~/.droptoken.rds")
	~~~

- Then pass the token to each drop_ function
	~~~
	d = drop_acc(dtoken = token)
	dropdir = "HBC Team Folder (1)/Consults/firstname_lastname/	$technology_of_$intervention_on_$tissue_in_$organism_$hbccode
	~~~

- All reports
	~~~
	sync(".", remote = dropdir, token = token, pattern = ".html")
	sync(".", remote = dropdir, token = token, blacklist = c(“Rproj”, “rda”, ...))
	~~~

2) By hand, by zipping up results and copying them to the appropriate folder on Dropbox.
It’s a good idea to always include the code you used for the results, as well as any linked results
[*FUTURE: NEED DISCUSSION*]

**We don’t use Dropbox to share bcbio results, BAM files, fastqs, or bcbio objects. If people needs those, we point them to the server or have them come with harddrive.**


## Suggested structure for the project folder that is in git repository or dropbox
~~~
Analysis
config
metadata
docs
templates
README
reports (code goes to GIT REPO, DROPBOX IF YOU WANT)
RMD (go to DROPBOX)
HTML (go to DROPBOX)
Data (R objects that you don’t want to sync to any place)
Results (go to DROPBOX) dropn[#FUTURE: NEED DISCUSSION]
~~~




<!--stackedit_data:
eyJoaXN0b3J5IjpbLTE0MzIxMDA1OTFdfQ==
-->
