# Git tips

- [Pro git book](https://git-scm.com/book/en/v2)
- https://github.com/Kunena/Kunena-Forum/wiki/Create-a-new-branch-with-git-and-manage-branches
- https://nvie.com/posts/a-successful-git-branching-model/
- http://sandofsky.com/blog/git-workflow.html
- https://blog.izs.me/2012/12/git-rebase
- https://benmarshall.me/git-rebase/
- [find big files in history](https://stackoverflow.com/questions/10622179/how-to-find-identify-large-commits-in-git-history)
- [remove a big file from history](https://www.czettner.com/2015/07/16/deleting-big-files-from-git-history.html)


# merge master branch into (empty) main and delete master
```
module load git
git fetch origin main
git branch -a
git checkout main
git merge master --allow-unrelated-histories
git add -A .
git commit
git push

git branch -d master
git push origin :master
```

# Add remote upstream
```
git remote -v
git remote add upstream https://github.com/ORIGINAL_OWNER/ORIGINAL_REPOSITORY.git
```
# Sync with upstream/master, delete all commits in origin/master
```
git fetch upstream
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
# big feature workflow - rebase - squash
```
# sync master with upstream first
# create new branch and switch to it
git checkout -b feature1
# create many commits with meaningfull messages
git add -A.
git commit
git push
# upstream accumulated some commits
git fetch upsteam
# rebasing the branch not the master
# to PR from the branch later not from the master
git rebase upstream/master
# see latest commits from HEAD down to the start of feature1
# on top of upstream
git log --oneline --decorate --all --graph
# interactive rebase for the last 13 commits (including head)
git rebase -i HEAD~13
# set s (squash) in the interactive editor for all commits except for the top one
# alter commit message
# force since origin as 13 separate commits
git push --force origin feature1
# PR from feature1 branch to upstream/master
```

# 2 Feature workflow
```
git checkout -b feature1
git add -A.
git commit
git push --set-upstream origin feature1
# pull request1
git checkout master
git checkout -b feature2
git add -A.
git commit
git push --set-upstream origin feature2
# pull request 2
```

# Feature workflow w squash
```
git checkout -b feature_branch
# 1 .. N
git add -A .
git commit -m "sync"

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

# get commits from maintainers in a pull request and push back
```
git fetch upstream pull/[PR_Number]/head:new_branch
git checkout new_branch
git add 
git commit
git push --set-upstream origin new_branch
```

# ~/.ssh/config
```
Host github.com
    HostName github.com
    PreferredAuthentications publickey
    IdentityFIle ~/.ssh/id_rsa_git
    User git
```

# Migrating github.com repos to [code.harvard.edu](https://code.harvard.edu/)

See [this page](https://gist.github.com/niksumeiko/8972566) for good general guidance

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
# commit saved copy of download_data.md without secrets
```
