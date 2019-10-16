# Git tips

- [Pro git book](https://git-scm.com/book/en/v2)
- https://github.com/Kunena/Kunena-Forum/wiki/Create-a-new-branch-with-git-and-manage-branches
- https://nvie.com/posts/a-successful-git-branching-model/
- http://sandofsky.com/blog/git-workflow.html
- https://blog.izs.me/2012/12/git-rebase

# Sync with upstream/master, delete all commits in origin/master branch
```
git checkout master
git reset --hard upstream/master
git push --force
```

# Sync with upstream/master
```
git fetch upstream
git checkout master
git merge upstream/master
```

# Feature workflow
```
git checkout -b feature_branch
# 1 .. N
git add -A .
git commit -m "sync"
git push?

git checkout master
git merge --squash private_feature_branch
git commit -v
git push
# pull request to upstream
# code review
# request merged
git branch -d feature_branch
git push origin :feature_branch
```

# Migrating github.com repos to [code.harvard.edu](https://code.harvard.edu/)

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

# Remove sensitive information from the file and from the history
```
Make a backup
# cd ~/backup
# git clone git@github.com:hbc/knowledgebase.git
cd ~/work
git clone git@github.com:hbc/knowledgebase.git
git filter-branch --tree-filter 'rm -f admin/download_data.md' HEAD
git push --force-with-lease origin master
# commit a saved opy of download_data.md without secrets
```
