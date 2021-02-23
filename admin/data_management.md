# Data managment on O2

As an HMS Core we get storage on O2 for free. We are not the biggest user but we are in the top 10. As such it is smart for us to be good citizens and keep our footprint as low as possible. Ways to reduce our footprint include
## Reduce space used by active consults
* delete redundant files 
  * a special case of this is when you receive raw data as some seq facilities (cough, BPF, cough) will send us both foo.fastq.gz and foo.fastq.bz2 files, we can delete the bz2 files
* run bcbio analyses so temporary work folders are put on scratch 
  * the brute force approach to this is to setup on scratch with everything as a symlink, run the analysis and move the final folder over to the PIs project folder
  * a more elegant way to do this is to set the output folder of the bcbio run to be in the project folder in the PIs folder
* keep project folders tidy, delete things you are no longer using (this is a judgement call and not really enforced but can be an issue once the project is complete)
* compress raw data files
  * bcbio can handle gzipped and bzipped (bz2) files, compress those raw fastqs!

## Get the data out of our main storage area
### Return data from completed analyses to the researcher
* this is a hassle so we typically only do it for larger data folders
* can be by
  1) "sneaker net" - downloading onto a drive and handing it off to the researcher (not preferred)
  2) GLOBUS - https://www.globus.org/ If the researcher sets up a GLOBUIS account, we can work with HMS RC to make their directory available to GLOBUS so the researcher can download it through the web or CL interface. While we haven't tested it, this method has the potential to be the least work for us.
  3) upload to the researcher's server - preferred if access to the server is simple. A good example would be the research data storage that HMS PIs have access to or an FTP server. Passwords, logins and occasionally VPN access are typically required.
  
I recommend avoiding things like Dropbox, Google Drive or Box unless the data is small. They aren't really built for this purpose.

### Archive the data
We have access to standby storage on O2 (/n/standby/cores/bcbio/). For projects that are either too small to bother with returning to the researcher or projects where we think we may want to access the data again, we can tar.gz them and store them here. Leave a symlink in the original directory to allow easy restoration of the project 

**Once you have restored the project, delete the standby file.**

**Once you are finished with the project, rearchive it**

Please don't keep an archived copy of the diretory in two places plus the expanded folder. Duplicated data is wasted space and makes John cry.

## How to decide what to do with the data
With the caveat that every project is different here are some general guidelines to help guide your decision making process.

1) Is the data "large" (>500GB)? As much as possible, we'd prefer to get rid of these ASAP

2) Will you need to access the data and derived files again? If yes, tidy up any unnecessary files and archive
 
 The following points can inform your decision making about how likely we will need to reaccess the data
  * Is the analysis published? If its published, its likely we won't be using it again. Ask the researcher what they want returned to them and delete the rest.
  * If the data is unpublished, did the analysis work? If the data is garbage, we likely won't be using it again. Ask the researcher what they want returned to them and delete the rest.
  * Did the consult end well? If it didn't its likely we won't be using it again. Ask the researcher what they want returned to them (with a time warning in case they don't respond) and delete the rest.
  * Is the consult unique in its approach? i.e. do you think its something we might come back to in the future for other analyses. If yes, tidy it up and archive. If it is something like a standard RNA-seq analysis, we probably don't need it. 
 * How long has it been since the data was accessed? If its older than 2 years, we likely won't be using it again. Ask the researcher what they want returned to them and delete the rest.

## Globus howto (in progress)

1. ```tar cf project_dir.tar project_dir```
2. login to www.globus.org
3. File Manager / Collection / type HMS-RC
