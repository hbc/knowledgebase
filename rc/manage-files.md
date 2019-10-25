---
title: Managing files
description: This code helps with managing file names
category: computing
subcategory: tips_tricks
tags: [bash, osx, linux]
---

## How to remove all files except the ones you want:

First, expand rm capabilities by:
`shopt -s extglob`

Then use find and remove:
`find . ! -name 'file.txt' -type f -exec rm -f {} +`


## Rename files

  It has a lot of good options but one pretty useful for removing whitespaces and set all to lowercase:

  `rename -c --nows <fileName>`
 
## Use umask to restrict default permissions

Set umask 007 in our .bashrc. Then newly created directories will have 770 (rwxrwx---) permissions,
and files will have 660 (rw-rw-000).
