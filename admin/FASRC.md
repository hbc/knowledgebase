# FASRC

Login:

From terminal window:

ssh userid@login.rc.fas.harvard.edu

Enter pw (you won’t see the cursor move)
Enter verification code from DUO (mine is labeled “Third-Party msimoneau@login.rc.fas.harvard.edu”) (again, cursor won’t move)

Note - your userid is typically the first initial of your first name and your last name. (ie. jsmith)

You’ll land in your home directory (ie /n/home04/userid). Your home directory provides 100GB of storage that is private to you. It cannot be expanded or shared. It is not suitable for storing data security Level 3 or above. 

The paths for shared work/storage are: 

/n/holylfs05/LABS/hsph_bioinfo or 
/n/netscratch/hsph_bioinfo

We are in the group hsph_bioinfo. Three folders are available: /n/holylfs05/LABS/hsph_bioinfo (Everyone, Lab, Users directories). The hsph_bioinfo group has a 20T storage limit. Additionally,, /n/netscratch/hsph_bioinfo has 50TB total available for the group. Files are not backed up and will be removed after 90 days. Use this storage if you’re running low elsewhere.

Primarily everyone should use the Lab Directory. (/n/holylfs05/LABS/hsph_bioinfo/Lab and /n/netscratch/hsph_bioinfo)

drwxrwsr-x+  5 root hsph_bioinfo 4096 Mar  6 11:44 Everyone

drwxrws---+  5 root hsph_bioinfo 4096 Aug  9  2021 Lab

drwxrws---+ 12 root hsph_bioinfo 4096 Apr  9 14:18 Users

- Everyone – This directory is set drwxrwsr-x and is useful for sharing data with Other (all other cluster users). This directory is not available from Globus. Bcbio was here in the past. Could be used for nextflow scripts, etc.
- Lab – Using the Lab directory is the preferred method to share data with your lab group. This directory is also visible to the owner (and only the owner) when accessing the share via Globus (if applicable).
- The Users directory only exists on shares created prior to October 2024. On new shares and /n/netscratch please use the Lab directory to create personal folders.  The sub-directories in the Users directory are owned by individuals

There’s an FASRC group named hsph_bioinfo_admins - Shannan and Lorena are in this group and have access to everything in /n/holylfs05/LABS/hsph_bioinfo

## Copy data to/from FASRC using scp, ftp, Globus, etc. 

**SCP**

Zhu: i use scp to copy files from O2. 

Cd to destination folder on FAS then:

scp -r userid@transfer.rc.hms.harvard.edu:path/to/your/folder .


**Globus**

- use Harvard FAS RC Holyoke for access to /n/holylfs05/LABS/hsph_bioinfo/ and netscratch - 2025
  
/n/holylfs05/LABS/hsph_bioinfo - Just Lab and Users directories are available from globus. Copy data to these two directories (or better yet netscratch.) The “Everyone” directory is to store data available to everyone on FASRC, but is not exposed to the wider Globus world.

Move data for transfer from Everyone to Lab, then it will be visible in Globus.

lftp - for GEO upload

lftp is one of the options GEO allows for data upload Use the -c option to make the upload continue if the connection drops for some reason

Adam/FASRC was helpful for a GEO upload Oct 2024 

- use lftp  -c to keep an ftp running!!! 

do the ftp from a FASRC login node

## Misc Notes:

- holyscratch01 will be set to read-only during this december maintenance and will be decommissioned February 1, 2025. Please move any needed scratch data to netscratch and begin using it instead if you have not done so already. The global $SCRATCH variable will be changed to /n/netscratch
- More on FASRC storage: https://docs.rc.fas.harvard.edu/kb/quickstart-guide/  -- Determine where your files will be stored
- FASRC - use chown command to change ownership of files. Ownership doesn’t matter for billing purposes.
