# Guidelines on downloading from different sequencing facilities

## Biopolymers

bipolymers make their data available through an FTP site.

Their data will often come as both fastq files AND fastq.bz2 files. Feel free to delete one of these, we really don't need both and they take up a lot of space.


## Dana Farber MBCF

### Notes
- Zach and co. always share raw data but may default to sharing it through their pydio web interface, which is not reliable. 

- If you email Zach (zherbert@mail.dfci.harvard.edu) or Drew (andrew_caruso@mail.dfci.harvard.edu) and tell them who's data you need (cc: the researcher) they will setup an FTP site for you to use. 

- Make sure to let them know you've pulled down the data, so they can turn off the site when you're done (it costs money to run this).

- Their data is typically in tar.gz files, it can pay off to decompress them right away so you know if you have the whole file...

### Getting the data

- You can access the data through a `wget` command. 
- Preface it with `nohup` so your job keeps running even if you connection drops. 
- The final nohup.out file will the download progress in it if you want to confirm download. 
A typical command might be something like this:

`nohup wget -m ftp://jhutchinson:MBCFjhutchinson\!@34.198.31.178/*`

*note the escaped exclamation point in the password, they like to put characters like that in their passwords*

## Broad Institute

It can depend on the platform the researcher used, but the Broad typically only give out BAM files for normal RNA-seq runs. For their DGE (96 well) platform, they give everything under the sun.
You have to use their ASPERA system to pull down the files and you will need not only login and password info to get the data, but a limited time password to decrypt the data. While you can do this on your machine but by far the easiest method [*can someone update this when they do an actual Broad download? I'm not sure this is toally accurate*] is with their command line script `shares_download.sh`. If you don't have this info, have the client request it from the Broad and have the client send you the Broad's reply email with the info.


