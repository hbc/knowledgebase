---
title: Data downloading
description: Guidelines on downloading from different sequencing facilities.
category: admin
subcategory: guides
tags: [tutorial]
---

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

- If you email Zach (zherbert@mail.dfci.harvard.edu) or Drew (andrew_caruso@mail.dfci.harvard.edu) and tell them who's data you need (cc: the researcher) they will setup an FTP site for you to use.

- Make sure to let them know you've pulled down the data, so they can turn off the site when you're done (it costs money to run this).

- Their data is typically in tar.gz files, it can pay off to decompress them right away so you know if you have the whole file...

## Getting the data

- You can access the data through a `wget` command.
- Preface it with `nohup` so your job keeps running even if you connection drops.
- The final nohup.out file will the download progress in it if you want to confirm download.
A typical command might be something like this:

`nohup wget -m ftp://jhutchinson:MBCFjhutchinson\!@34.198.31.178/*`

*note the escaped exclamation point in the password, they like to put characters like that in their passwords, which are usually in the form of MBCF$userid!*

# Broad Institute

- It can depend on the platform the researcher used, but the Broad typically only give out BAM files for normal RNA-seq runs.
- For their DGE (96 well) platform, they give everything under the sun.
- You have to use their ASPERA system to pull down the files and you will need not only login and password info to get the data, but a limited time password to decrypt the data.
  - you can run ASPERA on your machine but by far the easiest method [*can someone update this when they do an actual Broad download? I'm not sure this is toally accurate*] is with their command line script `shares_download.sh`.
  - If you don't have this info, have the client request it from the Broad and have the client send you the Broad's reply email with the info.

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

# Basespace by Sergey
basespace-cli. New Illumina sequencers upload data to basespace cloud. bs utility copies data from cloud to HPC. To copy bcl files: bs cp //./Runs/<project_name>/Data .
