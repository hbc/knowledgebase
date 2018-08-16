
# Contributing to knowledegebase

You can contribute to knowledgebase in one of two ways. 

1. Add the information as **an issue**
2. Add the information via a **pull request**

## Add the information as **an issue**

* 
* 

## Add the information via a **pull request**

* `git clone <repo>`
* `git checkout master`
* `git pull origin master`
* `git checkout -b <branch_name>`
* *make necessary changes locally*
* `git add <changes>`
* `git commit -m "good message"`
* `git push origin <branch name>`
* *create pull request in github to merge the new branch with master branch*

After the pull request is created, one other member with write edits will have to "review" it, at which time they can approve it first, then merge it! 
**The person who merges should make sure they *delete* the branch they pulled from!**

In both cases, you will be required to **tag what you are adding with the appropriate labels**. For every issue/ pull request you can add tags based on the following questions:

### Is this information you would like to add to the knowledgebase or is it an internal feature request for the develeopers?

Tag options (must choose one): `docs` or `feature request`


### Which major category (Level 1) does this fall under?

Tag options (must choose one): `admin`, `research`, `computing`, `faq`

### Which sub-category does this classify with?

Each Level 1 category has the following Level 2 tags (must choose one):

* Admin: `guides`, `tools`, `internal resources`, `general resources`
* Computing: `bcbio`, `installation`, `tips and tricks`
* Research: `rnaseq`, `chipseq`, `atacseq`, `scrnaseq`, `smallrnaseq`, `general_ngs`, `wes`, `wgs`, `orphans`

#### Tag your issue with any one or more of the relevant tags

* `hpc`, `bcbio`, `local`, `R`, `python`, `snakemake`, `tutorial`, `template`, `bash` (more to come?)

