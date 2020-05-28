Main guide is here:
https://www.ncbi.nlm.nih.gov/geo/info/submission.html
Highh throughput sequencing is here:
https://www.ncbi.nlm.nih.gov/geo/info/seq.html

# The preparation 

## Analyst responsibilities
You will need
1) [GEO metadata sheet](https://www.ncbi.nlm.nih.gov/geo/info/seq.html)
2) Raw fastq files
3) Derived files for data  
  a) RNAseq
- raw counts table (as tsv/csv), can put as supplementary file
- TPM (as tsv/csv), can put as supplementary file
- bams are OK too but I have never been asked for them

4) details on the analysis for the GEO metadata sheet, including
- which sequencer was used
- paired or single end reads?
- insert size if paired
- programs and versions used in the analysis, including the bcbio and R portions

Example metadata sheets can be found in this Dropbox folder:
https://www.dropbox.com/sh/88035zd8h9qhvzh/AACmHB7xsXhdgrSyZY42uwLYa?dl=0

For all of the raw and derived data files, you will need to run md5 checksums.

## Researcher responsibilities
The client wil need to give you the details about things that were involved in the experiment and library preparation.
These include
- growth protocol
- treatment protocol 
- extract protocol
- library construction protocol
They will also need to supply the general info about the experiment including:
- title
- summary
- overall design
- who they want to be a contributor

I usually fill out what I can and then send them the metadata sheet with their areas to fill out highlighted.

# The upload
Once you have the data, derived data and metadata sheet, its time to upload to the GEO FTP server.
Sign into your NCBI and GEO account and go to the [Transfer Files](https://www.ncbi.nlm.nih.gov/geo/info/submissionftp.html) link on the GEO submission page. 

There they will tell you what your directory is on the GEO FTP server (for example, uploads/jnhutchinson_AtsZaoGM) as well as the server address (eg. ftp-private.ncbi.nlm.nih.gov) login (geoftp) and password (rebUzyi1). 

Go to your  upload directory on O2 with the GEO submission files and login to the ftp server using lftp geoftp:rebUzyi1@ ftp-private.ncbi.nlm.nih.gov: . Note that lftp is not available on login or interactive nodes, so you will need to ssh to the O2 transfer node (ssh user@transfer.rc.hms.harvard.edu) to use it. *You can also use Filezilla if your files are on your local machine.* Then move to your remote upload directory (cd /uploads/jnhutchinson_AtsZaoGM, *for Filezilla, you should enter this into the Remote site: directory box*) and start your upload. For lftp, you can use 
```mirror -R``` or ```mput *``` to upload the files. For Filezilla, just drag the files over to the remote directory. The sit back and maybe work on something else, or like, take a break from the bionformatics mine while everything uploads. If you have a ton of files, you may want to use something like tmux to prevent your session from being terminated. 


When the upload is complete, notify GEO of the submission using the cleverly named [Notify GEO](https://submit.ncbi.nlm.nih.gov/geo/submission/) link. 

You will receive an email confirming your upload and GEO staff will contact you if there are any issues. Common issues to watch out for are:
1) column headings in derived data not matching fastq sample names
2) missing gene ids in derived data

Less commonly they may ask you to fix insufficiently descriptive summary and overall design

# The aftermath

Note that unless you specifically set things up otherwise, the submission will be tied to your name and you will have to be responsible for updates and releases (i.e. you will be the "Investigator"). You can deal with this one of two ways
1) set things up from your initial login to have you as the submitter and the researcher as the Investigator, I personallly find this inconvenient as I may be doing multiple GEO submission for different researchers, but YMMV
2) do the submisison yourself as the Investigator and submitter and once the submisison is accepte, email GEO to have the submission transferred to the researcher. 
Note that both of these methods will require the researcher to obtain both an NCBI (if they don't already have one) account and share the login id and email address with you. 
