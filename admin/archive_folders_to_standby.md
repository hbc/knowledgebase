# archive folders

## 1) Login to specified login node on O2, via ssh
In the example below, my login ID is `jnh7` and the login node I am using is `login05.o2.rc.hms.harvard.edu` 
You will need to change `jnh7` to your ID to login. 

Example command:
`ssh -XY -l jnh7 login05.o2.rc.hms.harvard.edu`

Login with your usual password and DUO challenge

## 2) Run tmux
Running tmux will let you disconnect from your “session” while still letting your command run in the background.. This is super useful if the command is going to take a long time to run and you need to shut down your computer or disconnect form the wifi for some reason

The command below will setup a named tmux session called 	`foo`. Change `foo` to something that will help you remember what you are working on in this tmux session.

Example command:
`tmux new -s foo`


## 3) Start an interactive session

You should run the folder compression in an interactive session instead of the login node as the login nodes are shared and running commands in them can slow down the cluster for everyone. As such, RC may kill commands that take too much memory on login nodes without notifying you.

The command below will give you an interactive node with 8 gigs of RAM (`-mem 8000M`) for 8 hours (`-t 0-08:00`)

Example command:
`srun --pty -p interactive --mem 8000M --x11 -t 0-08:00 /bin/bash`


## 4) Compress the folder
Go to the folder on O2 
In this example I am pretending to work on  the project we did with Arlene Sharpe under the hbc code hbc03895 
This project is in the following folder:
`/n/data1/cores/bcbio/PIs/arlene_sharpe/Sharpe_RNAseq_analysis_of_siRNA_treated_PANCafs_and_myeloid_cells_after_coculture_hbc03895`

To compress it I can issue the following  commands
	Go to the PI  folder:
`cd /n/data1/cores/bcbio/PIs/arlene_sharpe/`
	Compress the project folder with tar and gzip in a single command and 
a) use the same project name as part of the compressed file name **AND**
b) add the date to the compressed file (I’ll use the date I wrote this document in the format YYYYMMDD, or 20201020)
c) to ensure the compression finished,  add the `—remove-files` option to the tar command, this will remove the files once the compression has scceeded. 
In general, the command to compress files would then look like this (`--remove-files` has to be before the other options as `-f` indicates a file name coming after):
`tar --remove-files -cvzf  YYYYMMDD_folder.tar.gz folder`

And in the case of this Arlelen Sharpe hbc038956  example, it would look like this:
`tar --remove-files -cvzf  20201020_Sharpe_RNAseq_analysis_of_siRNA_treated_PANCafs_and_myeloid_cells_after_coculture_hbc03895.tar.gz Sharpe_RNAseq_analysis_of_siRNA_treated_PANCafs_and_myeloid_cells_after_coculture_hbc03895`
This will compress the folder  and remove the original files (which will still be around in the backed up .snapshot folder for some weeks or months!)

### Additional considerations:


The `—remove-files` option ... The files are not *supposed* to be removed unless the tar operation is successful, however this apparently is not completely bulletproof!  An alternative suggestion is to not use -remove-files in the tar command. After checking that everything is correct (see subsequent suggestions), then delete the original files.

Timing out. It is possible, for some very large directories, for an srun of even 12 hours to time out before compression is complete. If you don't notice this, or if the complete compression operation fails for some other reason, there *will* be an apparent tar.gz file created, though incomplete or corrupted.  You could easily copy this over to standby and never realize the file was corrupt. So you should check for file integrity before copying.

`gzip -t file.tar.gz` This uses gzip's built‑in test mode to verify the compressed data's checksum. If the file is corrupted or incomplete, gzip will report an error.

Examine checksums of the compressed file before copying to standby (see 6 below) and the compressed file after copying to standby.  If the checksums match, you can safely delete the original data and original compressed file.

`md5sum file.tar.gz` will return the checksum.

Consider removing large raw data files, like .bam and .fastq files before compressing the directory. 

List them first. For example, `find . -type f -name '*.bam'` will list the .bam files in the directory and subdirectories.

But *only* delete these if you have written confirmation from the client that they have copies of these (and if you trust the client to have copies of these.) It is most likely that you are archiving a project after it has been completed, and maybe a GEO submission has already been made. But if a GEO submission is forthcoming, even if the client says they have copies of the raw data, I might wait and not delete these until after the GEO submission has been completed. 

`find . -type f -name '*.bam' -delete` will delete them.


## 5) Login to the transfer node
The standby folder can only be accessed from the transfer node

SSH to the transfer node
`ssh transfer`
Login with your password and DUO challenge, it shouldn’t require your login ID

## 6) Move  compressed folder to the standby folder

Our standby folder is located at `/n/standby/cores/bcbio/compute/archived_reports/tier2`
You can see that there are already a bunch of tar gzipped folders in there. 

Move your newly compressed folder to the standby folder using rsync.
Srync follows the general pattern of
`rsync -options source destination`

I typically use the following rsync options (as `-ravzuP`)
* r = recursive, ie. transfer all the subfolders too
* a = archive, preserve everything
* v= verbose, print  updates during tnansfer to screen
* u = update, skip files on receiver that are newer (this is useful if we screwup and have already archived the folder)
* P = show estimate progress as a bar or percentage transferred

Here I will continue with the example from 5 using Arlene Sharpe’s analysis using our compressed folder as source and the standby folder as destination

Example command:
`rsync -ravzuP /n/data1/cores/bcbio/PIs/arlene_sharpe/20201020_Sharpe_RNAseq_analysis_of_siRNA_treated_PANCafs_and_myeloid_cells_after_coculture_hbc03895.tar.gz /n/standby/cores/bcbio/compute/archived_reports/tier2/`

This will copy/sync the compressed folder over to the standby folder. 

Once the file is done copying/syncing, you can erase the source file

Example command:
`rm /n/data1/cores/bcbio/PIs/arlene_sharpe/20201020_Sharpe_RNAseq_analysis_of_siRNA_treated_PANCafs_and_myeloid_cells_after_coculture_hbc03895.tar.gz`

## 7) Link the standby copy of the compressed folder back to the PI folder
This should help locate it again in the future
We do this using a “symlink”, which is similar to a “Shortcut” in Windows or and “Alias” in OSX. 
The command to make a symlink in Unix is `ln -s` and typically has the format of `ln -s  source destination`

Using our example of the Sharpe analysis, the command would be:

`ln -s /n/standby/cores/bcbio/compute/archived_reports/tier2/20201020_Sharpe_RNAseq_analysis_of_siRNA_treated_PANCafs_and_myeloid_cells_after_coculture_hbc03895.tar.gz /n/data1/cores/bcbio/PIS/arlene_sharpe/ 20201020_Sharpe_RNAseq_analysis_of_siRNA_treated_PANCafs_and_myeloid_cells_after_coculture_hbc03895.tar.gz`

This will drop a symlink in the arelene_sharpe PI directory. W
You can see that you have a symlink by running `ls -lh ` in the PI directory, you should see at least one line with  the `20201020_Sharpe_RNAseq_analysis_of_siRNA_treated_PANCafs_and_myeloid_cells_after_coculture_hbc03895.tar.gz` file listed with an arrow beside it pointing to the standby folder
i.e. `20201020_Sharpe_RNAseq_analysis_of_siRNA_treated_PANCafs_and_myeloid_cells_after_coculture_hbc03895.tar.gz -> /n/standby/cores/bcbio/compute/archived_reports/tier2/20201020_Sharpe_RNAseq_analysis_of_siRNA_treated_PANCafs_and_myeloid_cells_after_coculture_hbc03895.tar.gz`






## NOTES

If you get disconnected from your session , you will be able to go back to this session by logging in as above in step 1 again and running the command to reconnect to tmux session `foo`

Example command:
`tmux attach -t foo`
