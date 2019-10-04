---
title: Migrating legacy repos to code.harvard.edu
description: GitHub instructions and commands to migrate
category: admin
subcategory: guides
tag: [github, tutorial]
---


## Migrating old repos to Enterprise

These are instructions for individual repos you may have locally, for older projects that you would like to move over to Enterprise [code.harvard.edu](https://code.harvard.edu/)

1. Set up your ssh keys. You can use your old keys (if you remember your passphrase) by going to `Settings --> SSH and GPG keys --> New SSH key`
2. Create your repo in code.harvard.edu. Copy the 'Clone with SSH link`:  `git@code.harvard.edu:HSPH/repo_name.git` (*NOTE: some of us have had trouble with the HTTPS link*)
3. Go to your local repo that you would like to migrate. Enter the directory.

```
# this will add a second remote location
git remote add harvard git@code.harvard.edu:HSPH/repo_name.git

# this will get rid of the old origin remote
git push -u harvard --all
``` 

4. You should see the contents of your local repo in Enterprise. Now go to 'Settings' for the repo and 'Collaborators and Teams'. Here you will need to add Bioinformatics Core and give 'Admin' priveleges.


> **NOTE:** If you decide to compile all your old repos into one giant repo (i.e. [hbc_mistrm_reports_legacy](https://code.harvard.edu/HSPH/hbc_mistrm_reports_legacy)), make sure that you remove all `.git` folders from each of them before committing. Otherwise you will not be able to see the contents on each folder on Enterprise.
