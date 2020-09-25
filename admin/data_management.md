# Data managment on O2

As an HMS Core we get storage on O2 for free. We are not the biggest user but we are in the top 10. As such it is smart for us to be good citizens and keep our footprint as low as possible. Ways to reduce our footprint include
## Reduce space used by active consults
* delete redundant files 
  * a special case of this is when you receive raw data as some seq facilities (cough, BPF, cough) will send us both foo.fastq.gz and foo.fastq.bz2 files, we can delete the bz2 files
* run bcbio analyses so temporary work folders are put on scratch 
  * the brute force approach to this is to setup on scratch with everything as a symlink, run the analysis and move the final folder over to the PIs project folder
  * a more elegant way to do this is to set the output folder of the bcbio run to be in the project folder in the PIs folder
* keep project folders tidy, delete things you are no longer using (this is a judgement call and not really enforced but can be an issue once the project is complete)
\
## Get the data out of our main storage area
### Return data from completed analyses to the researcher
* this is a hassle so we typically only do it for larger data folders
* can be by
  1) "sneaker net" - downloading onto a drive and handing it off to the researcher (not preferred)
  2) GLOBUS - https://www.globus.org/ If the researcher sets up a GLOBUIS account, we can work with HMS RC to make their directory available to GLOBUS so the researcher can download it through the web or CL interface
  3) upload to the researcher's server - preferred if access to the server is simple. A good exmaple would be the research data storage that HMS PIs have access to. Passwords and logins are typically required.
I recommend avoiding things like Dropbox, Google Drive or Box unless the data is small. They aren't really built for this purpose.
### Archive the data
We have access to standby storage on O2 (/n/standby/cores/bcbio/). For projects that are either too small to bother with returning to the researcher or projects where we think we may want to access the data again, we can tar.gz them and store them here. Leave a symlink in the original directory to allow easy restoration of the project 
**Once you have restored the project, delete the standby file.**
**Once you are finished with the project, rearchive it**
Please don't keep an archived copy of the diretory in two places plus the expanded folder. Duplicated data is wasted space and makes John cry.
