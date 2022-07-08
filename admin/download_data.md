# Google drive - big file

Use rclone to copy or sync files to a Google Drive directory. Instructions are at https://rclone.org/drive/. A few notes:

- log in to the transfer node on HMS (`transfer.rc.hms.harvard.edu`)
- create new remote using `rclone config` and select `n`. Give it a name.
- choose `12` for Google Drive
- skip `client_id` and `client_secret`
- for `scope`, select either 1 if you want to modify on Google drive or 2 if you just want to read files and copy them
- no need to select a root unless you want access to the "Computers" folder
- skip `service_account_file`
- `n` for auto config since we are a remote machine
- copy the link into your browser. Log in to Google. Copy key and paste into terminal
- if the drive is in your "Shared" folder on Google Drive, configure as a Shared Drive.

Additional notes:

- type `rclone lsd remote:` to see the files on your Google drive directory
- if you want to see what in your "Shared with me" drive, use --drive-shared-with-me as in `rclone --drive-shared-with-me lsd remote:`
- note that "Shared" and "Shared with me" are two different drives.

To copy, use `rclone copy source:directory destination:directory`

## Alternate instructions
https://www.quora.com/How-do-I-download-a-very-large-file-from-Google-Drive

# ActiveMotif

- They will give us a long link like this: ftp://ftp.activemotif.com<https://urldefense.proofpoint.com/v2/url?u=ftp-3A__ftp.activemotif.com_&d=DwMF-g&c=WO-RGvefibhHBZq3fL85hQ&r=CiePDOg3jDpiyuOnMdRMAf55kb0y979UKdT9-l_8xx4&m=sG-xxxx_AhGbB1PxIiUTGzFWLjIwpK-ff9Oza6JZM&s=Nmf-t3MTzOiZQ0WcCpKsFd6PMprEvuMK-izroBVUwB0&e=>

- We just need to get in to "ftp://ftp.activemotif.com" (some browsers might not open this)

- log in with the user, password.

- click on the folders and get the links.

- get on the terminal, transfer node. 
```
wget --user=user --password=password {link_for_necessary_folders}
# 'wget -m' for a folder.
```

# Biopolymers

- Biopolymers (BPF) make their data available through an SFTP site

- Their data will often come as both fastq files AND fastq.bz2 files.
- Feel free to delete one of these, we really don't need both and they take up a lot of space.

- Typically, the researcher will have an email from BPF with their login, password and the dataset id (usually in the form of FC_$number)
- Complete instructions can be found [here](https://genome.med.harvard.edu/documents/illumina/BPFNGS_FTP_Instructions.pdf)
- You can use `scp` or `rsync` to pull down the files

A typical command might look like this:

`scp -r jmubel2@bpfngs.med.harvard.edu:./FC_03443 . `


# Dana Farber MBCF

## Notes
- Zach and co. always share raw data but may default to sharing it through their pydio web interface, which is not reliable.

- If you email Zach (zherbert@mail.dfci.harvard.edu) and tell him who's data you need (cc: the researcher) he will setup an FTP site for you to use. (Drew has left, so don't email him anymore.)

- Make sure to let them know you've pulled down the data, so they can turn off the site when you're done (it costs money to run this).

- Their data is typically in tar.gz files, it can pay off to decompress them right away so you know if you have the whole file...

## Getting the data

- You can access the data through a `wget` command.
- Preface it with `nohup` so your job keeps running even if you connection drops.
- The final nohup.out file will the download progress in it if you want to confirm download.
A typical command might be something like this:

`nohup wget -m ftp://userid:password\!@34.198.31.178/*`

-m mirror is to copy a mirror image of the directory/data including all files and subfolders

Use this if nohup isn't working. Double check the UN, PW and IP address as they change.:
`wget -m ftp://jhutchinson:MBCFjhutchinson\!@34.198.31.178`

*note the escaped exclamation point in the password, they like to put characters like that in their passwords, which are usually in the form of MBCF$userid!*

# Broad Institute

- It can depend on the platform the researcher used, but the Broad typically only give out BAM files for normal RNA-seq runs.
- For their DGE (96 well) platform, they give everything under the sun.
- For 10x they may give you cellranger count output data matrices, and may or may not include bam files or even fastq files, so check you have what you need.
- You have to use their ASPERA system to pull down the files and you will need not only login and password info to get the data, but a limited time password to decrypt the data.
  - you can run ASPERA on your machine but by far the easiest method [*can someone update this when they do an actual Broad download? I'm not sure this is toally accurate*] is with their command line script `shares_download.sh`.
  - If you don't have this info, have the client request it from the Broad and have the client send you the Broad's reply email with the info.

Getting the Aspera client installed to facilitate transfer directly to o2 using the command line is initially a bit of a chore, but once installed makes data transfers simple. 

Install the Aspera Connect client for Linux to a directory on o2.Transfer.RC.hms.harvard.edu as follows:

1.	Go here to the Aspera site [IBM Aspera Connect](https://www.ibm.com/aspera/connect/)

2.	Click on "see all installers"

3.	Check the box for a recent linux version

4.	Scroll down to Fix package location to get ftp credentials

5.	sftp with given userid and server location. For example: `sftp vRmfnWoc@delivery04.dhe.ibm.com`  
	It will prompt you for ftp password
	
6.	Download files with mget: `mget*`

	This will download:  
	- ibm-aspera-connect_4.1.1.73_linux.tar.gz
	- IBM_Aspera_Connect_4.1_User_Guide_for_Linux.pdf
	- IBM_ASP_CONNECT_V4.1.1_RN_EN.pdf

7.	Extract the tarball: `tar -xf ibm-aspera-connect_4.1.1.73_linux.tar.gz`  
	This will open a shell script ibm-aspera-connect_4.1.1.73_linux.sh

8.	Run the shell script: `bash ibm-aspera-connect_4.1.1.73_linux.sh`  
	This will install the Aspera Connect client into .aspera/connect under your home directory

9.	Add the executable to your path: `export PATH=~/.aspera/connect/bin:$PATH`
        
10.	Download the shares_download.sh script from http://www.broadinstitute.org/aspera/shares_download.txt:  
	`wget http://www.broadinstitute.org/aspera/shares_download.txt`
        
11.	Rename this bash script from shares_download.txt to shares_download.sh:  
	`mv shares_download.txt shares_download.sh`

12. Make it executible script: `chmod a+x shares_download.sh`
        
13.	Use the shares_download.sh script to download your file. You will need the shares site credentials from the email from the Broad with this usage:  
	`shares_download.sh /download/destination https://shares.broadinstitute.org SN0020420:password SN0020420/`
	
	For example, to download into a directory ./data :  
	`bash shares_download.sh data https://shares.broadinstitute.org SN0243649:OJCGRLNFB8O6R9P SN0243649/`
	
Steps 1:12 should only have to be run once.

# BaseSpace to O2 by Radhika (July 2022)

https://help.basespace.illumina.com/cmd-line-interfaces/basespace-cli/introduction-to-basemount#Overview

1. Log into the transfer node on O2: `ssh username@transfer.rc.hms.harvard.edu`
2. Use the command `basemount BaseSpace/` from wherever you want to mount (I mounted it in my home directory). If using it for the first time, you will have to authenticate using the link, see example below.
	```
	rsk27@transfer06:~$ basemount BaseSpace/
	,-----.                        ,--.   ,--.                         ,--.   
	|  |) /_  ,--,--. ,---.  ,---. |   `.'   | ,---. ,--.,--.,--,--, ,-'  '-. 
	|  .-.  \' ,-.  |(  .-' | .-. :|  |'.'|  || .-. ||  ||  ||      \'-.  .-'
	|  '--' /\ '-'  |.-'  `)\   --.|  |   |  |' '-' ''  ''  '|  ||  |  |  |  
	`------'  `--`--'`----'  `----'`--'   `--' `---'  `----' `--''--'  `--' 
	Illumina BaseMount v0.25.2.3271 public develop 2021-07-12 15:33
	
	Command called:
	    basemount BaseSpace/
	From:
	    /home/rsk27
	
	Mount point "BaseSpace/" doesn't exist
	Create this mount point directory? (Y/n) Y
	Creating directory "BaseSpace/"
	Starting authentication.
	
	You need to authenticate by opening this URL in a browser:
	  https://basespace.illumina.com/oauth/device?code=U7my2
	...
	It worked!
	Your identification has been saved.
	
	Mounting BaseSpace account.
	To unmount, run: basemount --unmount /home/rsk27/BaseSpace
	```
3. `ls BaseSpace/` will show you what is available to you
	```
	rsk27@transfer06:~$ ls BaseSpace/
	IAP  Projects  README  Runs  Trash
	```
5. Since it is mounted now, you can simply use cp or rsync, if you prefer, to copy over the necessary files/directories into the appropriate location. `cp ~/BaseSpace/Projects/BS_46-RNA_S-21-1766_GAP375/Samples/[A-Z]_*/Files/*gz .`
6. To unmount, run `basemount --unmount ~/BaseSpace`

# Basespace by Rory
wget https://da1s119xsxmu0.cloudfront.net/sites/knowledgebase/API/08052014/Script/BaseSpaceRunDownloader_v2.zip
unzip BaseSpaceRunDownloader_v2.zip
rm run_BaseSpaceRunDownloader.bat  BaseSpaceRunDownloader_v2.zip
python BaseSpaceRunDownloader_v2.py -r <Run ID> -a <access token>

if project is specified instead of run:
wget https://gist.githubusercontent.com/rlesca01/7ce2ca0c35c7ff97a215/raw/0eeaa8cc1b3eff00babf398a82a31f4b0946f5bb/BaseSpaceRunDownloader_v2a.py

# Basespace by Victor

Use Illumina's native GUI client or run [BaseMount](https://basemount.basespace.illumina.com) on Ubuntu.

The [Python downloader](https://support.basespace.illumina.com/knowledgebase/articles/403618-python-run-downloader) is deprecated and no longer supported by Illumina.

## Update:

If you need to download processed data (i.e: not the bcl files but the fastq or a whole project) or you don't have a run number you can use the [BaseSpace Sequence Hub CLI](https://developer.basespace.illumina.com/docs/content/documentation/cli/cli-overview#Entitiesandsub-entities). After installing it you can download specific datasets or projects with ```bs download```.

Example: Downloading specific fastq files from a project:

1st. Identify datasets with ```bs list datasets```
```
$bs list datasets > avail_datasets
$head avail_datasets

+-----------------------------+-------------------------------------+----------------------+---------------------+
|            Name             |                 Id                  |     Project.Name     |   DataSetType.Id    |
+-----------------------------+-------------------------------------+----------------------+---------------------+
| 13_Randa_01_25_19           | ds.9a309d19c84f44c191fc86919b9a562e | Randa_01_25_19_Run1  | common.files        |
| 11_Randa_01_25_19           | ds.c29329f59b2d42beaa0e617e50829e06 | Randa_01_25_19_Run1  | common.files        |
| 10_Randa_01_25_19           | ds.d3946986f8da4d8ba9a8c5db4037ff7c | Randa_01_25_19_Run1  | common.files        |
| 12_Randa_01_25_19           | ds.262b749568fb4733ad958f9dfb4df0e3 | Randa_01_25_19_Run1  | common.files        |
| 14_Randa_01_25_19           | ds.2645866a5f514e0fb0bb26b191eff138 | Randa_01_25_19_Run1  | common.files        |
| 23_Randa_01_25_19           | ds.e84152ac7b754e4ca4f0e859b5864344 | Randa_01_25_19_Run1  | common.files        |
| 5_Randa_01_25_19            | ds.3170cc50b2774650bdbcfedb5ce48830 | Randa_01_25_19_Run1  | common.files        |
```
(I clean the header and remove lines of ---- to make avail_datasets to dataset_list)

2nd. Clean up a bit the file and select the ones you're instered in. (in this case I'm using `gawk` to select the ones processed by "illumina.fastq.v1.8" analysis.
```gawk 'BEGIN{FS="|"}{print $2,$4,$5}' dataset_list | gawk 'BEGIN{OFS=","}$3=="illumina.fastq.v1.8"{print $1,$2}' > fastqIDs.txt```

Then I use the following code to download each dataset and store it in a Project folder:

`<obtain_ds.sh>`

```
#!/bin/sh

fastqID_file="$1"

while read -r line; do
    dsID="$(cut -d',' -f1 <<<"$line")"
    projectID="$(cut -d',' -f2 <<<"$line")"
	# line is available for processing
	bs download dataset -n ${dsID} -o ${projectID}/${dsID}
done < ${fastqID_file}
```

`<obtain_ds_job.sh>`
```
#!/bin/bash

#SBATCH --job-name=mittendorf          # Job name
#SBATCH --partition=short             # Partition name
#SBATCH --time=0-11:59                 # Runtime in D-HH:MM format
#SBATCH --nodes=1                      # Number of nodes (keep at 1)
#SBATCH --ntasks=1                     # Number of tasks per node (keep at 1)
#SBATCH --cpus-per-task=1             # CPU cores requested per task (change for threaded jobs)
#SBATCH --mem=4G                     # Memory needed per node (total)
#SBATCH --error=jobid_%j.err           # File to which STDERR will be written, including job ID
#SBATCH --output=jobid_%j.out          # File to which STDOUT will be written, including job ID
#SBATCH --mail-type=ALL                # Type of email notification (BEGIN, END, FAIL, ALL)


bash obtain_ds.sh fastqIDs.txt
```


# Basespace by Sergey
basespace-cli. New Illumina sequencers upload data to basespace cloud. bs utility copies data from cloud to HPC. To copy bcl files: bs cp //./Runs/<project_name>/Data .
