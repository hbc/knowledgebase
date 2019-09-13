# Guidelines for managing HBC Research Data

## Motivation
To handle data in a manner that allows it to be FAIR, i.e.

* Findable: associated with a unique identifier
* Accessible: easily retrievable
* Interoperable: "use and speak the same language" via use of standardized vocabularies
* Reusable: adequately described to a new user, have clear information about data usage, and have a traceable "owner's manual" or provenance

Many of the Core's standard operating procedures are geared towards reproducibility/reusability. Please adhere to the following guidelines.

**1. Where to analyze your data**

### High performance computing
Most of the data processing using bcbio is performed using the HMSRC cluster, O2 in the `/n/data1/cores/bcbio/PIs` directory. The exception is for FAS and HSPH faculty, where the preferred location is on the FASRC cluster, Odyssey. The exact location on Odyssey is TBD [Contact John].

When setting up an analysis, refer to the [Setting up an analysis guidelines](https://github.com/hbc/knowledgebase/blob/master/admin/guides/setting_up_an_analysis_guidelines.md) for how to name directories and which folders to include. 

As a reminder, it's important to make sure that you store a copy of the yaml config, metadata csv file, and slurm script used to run the analysis along with the raw data so that someone else can access the project and rerun it if necessary. It is also very helpful if you keep a copy of the bcb object for downstream analysis with your data. You should also commit the config, metadata and slurm script to code.harvard.edu.

### Downstream statistical analyses
Whether you perform your data analysis on the cluster or locally, it's important that the Rmd files be committed to code.harvard.edu as well as 


### O2

* The Core's shared space is located at `/n/data1/cores/bcbio/PIs`. 
* Refer to the [Setting up an analysis guidelines](https://github.com/hbc/knowledgebase/blob/master/admin/guides/setting_up_an_analysis_guidelines.md) for how to name directories and which folders to include. 
* Store a copy of the yaml config, metadata csv file, and slurm script used to run the analysis along with the raw data so that someone else can access the project and rerun it if necessary. 
* It is also very helpful if you keep a copy of the bcbio object for downstream analysis with your data.

### code.harvard.edu

Commit the following:
* metadata csv
* yaml config
* slurm submission script
* Rmarkdowns for analysis
* any auxillary scripts required to analyze the data

### Dropbox

Use dropbox to store results with collaborators
* html files
* Rmarkdown files used to generate the html files
* text files and "small" processed data files (i.e. files included in the reports, such as normalized counts)
* documents: manuscripts, extra metadata, presentations

### Basecamp

Use Basecamp for discussions of project progress.
* Link to Dropbox for results
* You may post small files on basecamp to share with collaborators
