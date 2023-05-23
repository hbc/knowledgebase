# Guidelines for managing HBC Research Data

[Docker cheatsheet](https://dockerlabs.collabnix.com/docker/cheatsheet/)

## Motivation
To handle data in a manner that allows it to be FAIR, i.e.

* Findable: associated with a unique identifier
* Accessible: easily retrievable
* Interoperable: "use and speak the same language" via use of standardized vocabularies
* Reusable: adequately described to a new user, have clear information about data usage, and have a traceable "owner's manual" or provenance

As a rule of thumb it may help to think of anything that is *only* on your local machine as being not reproducible or reusable.  

Many of the Core's standard operating procedures are geared towards reproducibility/reusability. Please try to adhere to the following general guidelines.


### O2
In general, O2 is for the big stuff. Also, anything needed to reproduce the results on run on the server should be here.
* The Core's shared space is located at `/n/data1/cores/bcbio/`.
  *  We keep projects for researcher in the PIs folder
    * this folder has the following structure 
    * PIfirstname_PIlastname/project_folder
    * ideally the project folder would match the github repo name and look similar to the Trello/Harvest name
    * it should at least contain the hbc code for tracking
* Refer to the [Setting up an analysis guidelines](https://github.com/hbc/knowledgebase/blob/master/admin/setting_up_an_analysis_guidelines.md) for how to name directories and which folders to include in the project folder. 
* Store a copy of the yaml config, metadata csv file, and slurm script used to run the analysis along with the raw data so that someone else can access the project and rerun it if necessary. 
* It is also very helpful if you keep a copy of the bcbio object for downstream analysis with your data.


### HBC org on github
In general, code is for a continually updated, searchable record of your code, and will mainly be made up of the type of code you run locally. 
Commit the following:
* Rmarkdowns for analysis
* any auxillary scripts required to analyze the data and/or present it in finished form (i.e. data object conversion scripts, yaml files for bcbioRNAseq, etc.)

### Dropbox

Use dropbox to share results and code with collaborators (`HBC Team Folder (1)/Consults/piFirstName_piLastName/project_folder`).
* html files
* Rmarkdown files used to generate the html files
* text files and "small" processed data files (i.e. files included in the reports, such as normalized counts)
* documents: manuscripts, extra metadata, presentations
* anything else the collaborator may need to reproduce the R-based analysis *IF* they wish to reproduce it
* *RNA-seq* example: create folders for QC, DE and FA. These may be in the base project folder, or if the project is complex, in folders labeled by dates. Within the QC folder, store the QC Rmd and the html report. In the DE folder, store the DE Rmd file, the DE report, and the raw and normalized counts matrices. In the FA folder, store the FA Rmd, the html report and the DESeq2 results tables with gene symbols added.

### Basecamp

Use Basecamp for discussions of project progress.
* Link to Dropbox for results
* You may post small files on basecamp to share with collaborators

### Cases to discuss:

* Many times, the client sends us a presentation (ppt/pptx) or a paper to better understand the project or to provide us with some information that may end up in the metadata. Where should we store these files?   
*John - a "docs" directory on O2 works for this*  
* Similarly, should we store original metadata files or only the csv file that we end up using.  
*John - I generally do, as it can contain important information or be easier to run past the original researcher. I mark it as "original" so I can clearly differentiate it from the metadata we ended up using. If available, I save any code I used to modify the original metadata with the metadata.*  
* Reviews. Consults that consist on reviewing the client paper (i.e code). Where to keep all the provided documents (if we have to keep them).  
*John - I could see keeping them on Dropbox, as we would typically want to share any reviewer responses with the researcher.*  
* In the case we have downloaded data from multiple flowcells, should we keep the originals or the concatenated files. Note that errors can happen when "preparing samples" (concatenating the wrong files for example).  
*John - would be great if we could only keep the concatenated files, but only after confirming no lane effects and proper concatenation.*
