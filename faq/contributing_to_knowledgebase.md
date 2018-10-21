---
title: Contributing to knowledgebase
description: This document describes procedures to follow when adding information to knowledgeBase
category: faq
subcategory: 
tags: [tutorial, github]
---


# Contributing to knowledegebase

You can contribute to knowledgebase in one of two ways. 

### (1) Add the information via a **pull request**

* `git clone <repo>`
* `git checkout master`
* `git pull origin master`
* `git checkout -b <new_branch_name>`
* make necessary changes locally
* `git add <changes>`
* `git commit -m "good message"`
* `git push --set-upstream origin <new_branch name>`
* continue to make necessary changes if need be (on github or locally)
* if contributing a document (markdown), ***annotate it appropriately*** (see below)
* create pull request in github ***annotated with appropriate labels*** to merge the new branch with master branch
* continue to make necessary changes if need be (on github or locally)

After the pull request is created, someone will have to "review" and "approve" it first, then merge it! 

**The person who merges should make sure they *delete* the branch they pulled from!**

### (2) Add the information as **an issue**

* Create an issue with a question or documentation 
* Tag with the necessary labels (see below)
* If you need to alert a specific person or persons, assign them to the issue

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

Please add all of the appropriate labels to tag an issue or a pull request being submitted based on the following guidelines:

* Add the `docs` label (mandatory).
* You will then need a single "category-" label (mandatory).
* You will then need a single "subcategory-" label (mandatory).
* You can then have as many other labels as necessary.


### Controlled vocabulary for headers/labels

#### "Category":

Use these as they are listed below when creating a header for a markdown file being created for contribution. When the pull request (or an issue) is created, use the matching label to tag it (matching labels will have the `category-` prefix).

* `admin` (label: `category-admin`)
* `research` (label: `category-research`)
* `computing` (label: `category-computing`)
* `faq` (label: `category-faq`)

#### "Sub-category" for specific "Categories":

Use these as they are listed below when creating a header for a markdown file being created for contribution. When the pull request (or an issue) is created, use the matching label to tag it (matching labels will have the `subcategory-` prefix).

* **Admin**: `guides`, `tools`, `internal resources`, `general resources`
* **Computing**: `bcbio`, `installation`, `tips and tricks`
* **Research**: `rnaseq`, `chipseq`, `atacseq`, `scrnaseq`, `smallrnaseq`, `general_ngs`, `wes`, `wgs`, `orphans`
* **Training**: `admin`, `materials`, `resources`

#### "Tags":

You can use as many of these as needed to populate the `tag:` item within the header of any markdown you are contributing. These are all also available as labels on github to tag pull requests and issues.

> *Note that these have been roughly split up into general, research and training so it is easier to read, but they are all interchangeable and can be used as needed.*

* **General tags**: `hpc`, `bcbio`, `local`, `R`, `python`, `snakemake`, `tutorial`, `template`, `bash`, `ssh`, `linux`, `osx`, `perl`, `report`
* **Research tags**: `quality_control`, `annotation`, `clustering`, `functional_analysis`, `differential_analysis`, `intron_retention`, `motif_analysis`, `trajectory`, `isoforms`, `visualization`, `readme`, `alignment`, `variant_calling`, `quantification`, `peak_calling`, `filtering`
* **Training tags**: `metrics`, `evaluation`, `catalyst`, `hsci`, `hms`, `hsph`, `social`, `literature`

#### Proposing a new label/"tag"

If none of the labels we have above are appropriate for your issue/pull request, you can propose a new label on #knowledgebase channel in the HBC team Slack.

* Tag file/issue/pull-request with the label you are proposing to add to the list.
* On the #knowledgebase channel in Slack, post a description of the (type of) file/issue/pull-request you are proposing it for, along with a 2 sentence justification for adding it. Please use "@channel" when posting this.
* The rest of the members of the core then will vote on it with a "thumbs up" or "thumbs down" reaction.
* There can be 3 outcomes from this:
	* Approved! - the post gets 3 thumbs ups (net value, so it could be 5 "up" and 2 "down")
	* Not approved! - everyone votes and final net value is 3 thumbs down
	* Discuss in group meeting - total net value is < 3 for either "up" or "down"
* If approved, move forward with adding the proposed label to master

