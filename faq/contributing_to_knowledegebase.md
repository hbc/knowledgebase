
# Contributing to knowledegebase

You can contribute to knowledgebase in one of two ways. 

1. Add the information as **an issue**
2. Add the information via a **pull request**

### Add the information as **an issue**

* 
* 

### Add the information via a **pull request**

* `git clone <repo>`
* `git checkout master`
* `git pull origin master`
* `git checkout -b <new_branch_name>`
* *make necessary changes locally*
* `git add <changes>`
* `git commit -m "good message"`
* `git push --set-upstream origin <new_branch name>`
* *create pull request in github to merge the new branch with master branch*

After the pull request is created, one other member with write edits will have to "review" it, at which time they can approve it first, then merge it! 
**The person who merges should make sure they *delete* the branch they pulled from!**

***

## Annotation

In both of the above cases, you will be required to **annotate what you are adding appropriately**. Below are a couple of ways to annotate your contribution; please note that they are not mutually exclusive:

### Annotate new markdown files before the pull request is created

* Add a yaml header with the 5 items listed in the following example:
```
---
title: Training Materials
description: Link to training materials
category: training
subcategory: resources
tag: [readme]
---
```

* "`title:`" and "`description:`" will obviously have the information specific to the file you are adding the header to. 
* For the remaining 3 items ("`category:`", "`subcategory:`", "`tag:`"), we have controlled vocabulary available and it is listed at the bottom of this section.
* Please note whatever you use as a "`tag:`" will need to be encased in brackets (`[]`).
* Make sure you only have one category and subcategory.
* However, you can use multiple tags separated by commas! E.g. `tag: [R], [visualization], [hpc]`

### Annotate issues and pull requests with existing labels

* Please add all of the appropriate labels to tag an issue or a pull request being submitted based on the following guidelines:
 ** Adding information to the knowledgebase: `docs` 
 ** Iternal feature request for the developers: `feature request`

### Which major category (Level 1) does this fall under?

Tag options (must choose one): `admin`, `research`, `computing`, `faq`

### Which sub-category does this classify with?

Each Level 1 category has the following Level 2 tags (must choose one):

* **Admin**: `guides`, `tools`, `internal resources`, `general resources`
* **Computing**: `bcbio`, `installation`, `tips and tricks`
* **Research**: `rnaseq`, `chipseq`, `atacseq`, `scrnaseq`, `smallrnaseq`, `general_ngs`, `wes`, `wgs`, `orphans`
* **Training**: `admin`, `materials`, `resources`

#### Tag your issue with any one or more of the relevant tags

* **General tags**: `hpc`, `bcbio`, `local`, `R`, `python`, `snakemake`, `tutorial`, `template`, `bash`, `ssh` 
* **Research tags**: `quality_control`, `annotation`, `clustering`, `functional_analysis`, `differential_analysis`, `intron_retention`, `motif_analysis`, `trajectory`, `isoforms`, `visualization`, `readme`, `alignment`, `variant_calling`, `quantification`, `peak_calling`, `filtering`

