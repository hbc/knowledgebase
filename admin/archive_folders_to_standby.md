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
In general, the command to compress files would then look like this:
`tar -cvzf --remove-files YYYYMMDD_folder.tar.gz folder`


And in the case of this Arlelen Sharpe hbc038956  example, it would look like this:
`tar -cvzf  20201020_Sharpe_RNAseq_analysis_of_siRNA_treated_PANCafs_and_myeloid_cells_after_coculture_hbc03895.tar.gz 20201020_Sharpe_RNAseq_analysis_of_siRNA_treated_PANCafs_and_myeloid_cells_after_coculture_hbc03895`
This will compress the folder  and remove the original files (which will still be around in the backed up .snapshot folder for some weeks or months!)

## 5) Login to the transfer node
The standby folder can only be accessed from the transfer node

SSH to the transfer node
`ssh transfer`
Login with your password and DUO challenge, it shouldn’t require your login ID

## 6) Move  compressed folder to the standby folder

Our standby folder is located at `/n/standby/cores/bcbio/archived_reports/tier2`
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
`rsync -ravzuP /n/data1/cores/bcbio/PIs/arlene_sharpe/20201020_Sharpe_RNAseq_analysis_of_siRNA_treated_PANCafs_and_myeloid_cells_after_coculture_hbc03895.tar.gz /n/standby/cores/bcbio/archived_reports/tier2`

This will copy/sync the compressed folder over to the standby folder. 

Once the file is done copying/syncing, you can erase the source file

Example command:
`rm /n/data1/cores/bcbio/PIs/arlene_sharpe/20201020_Sharpe_RNAseq_analysis_of_siRNA_treated_PANCafs_and_myeloid_cells_after_coculture_hbc03895.tar.gz`

## 7) Link the standby copy of the compressed folder back to the PI folder
This should help locate it again in the future
We do this using a “symlink”, which is similar to a “Shortcut” in Windows or and “Alias” in OSX. 
The command to make a symlink in Unix is `ln -s` and typically has the format of `ln -s  source destination`

Using our example of the Sharpe analysis, the command would be:

`ln -s /n/standby/cores/bcbio/archived_reports/tier2/20201020_Sharpe_RNAseq_analysis_of_siRNA_treated_PANCafs_and_myeloid_cells_after_coculture_hbc03895.tar.gz /n/data1/cores/bcbio/PIS/arlene_sharpe/ 20201020_Sharpe_RNAseq_analysis_of_siRNA_treated_PANCafs_and_myeloid_cells_after_coculture_hbc03895.tar.gz`

This will drop a symlink in the arelene_sharpe PI directory. W
You can see that you have a symlink by running `ls -lh ` in the PI directory, you should see at least one line with  the `20201020_Sharpe_RNAseq_analysis_of_siRNA_treated_PANCafs_and_myeloid_cells_after_coculture_hbc03895.tar.gz` file listed with an arrow beside it pointing to the standby folder
i.e. `20201020_Sharpe_RNAseq_analysis_of_siRNA_treated_PANCafs_and_myeloid_cells_after_coculture_hbc03895.tar.gz -> /n/standy/cores/bcbio/archived_reports/tier2/20201020_Sharpe_RNAseq_analysis_of_siRNA_treated_PANCafs_and_myeloid_cells_after_coculture_hbc03895.tar.gz`






## NOTES

If you get disconnected from your session , you will be able to go back to this session by logging in as above in step 1 again and running the command to reconnect to tmux session `foo`

Example command:
`tmux attach -t foo`