---
title: Tips for bcbio
description: This code helps with fixing jobs which timeout, and other general tips
category: computing
subcategory: bcbio
tags: [bcbio, bash, hpc]
---

## Installing a private bcbio development repository on O2
```bash
wget https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/scripts/bcbio_nextgen_install.py
python bcbio_nextgen_install.py ${HOME}/local/share/bcbio --tooldir=${HOME}/local --nodata
ln -s /n/app/bcbio/biodata/genomes/ ${HOME}/local/share/genomes
ln -s /n/app/bcbio/biodata/galaxy/tool-data ${HOME}/local/share/galaxy/tool-data
```

## How to fix potential conda errors during installation
Add the following to your `${HOME}/.condarc`:
```yaml
channels:
  - bioconda
  - defaults
  - conda-forge
safety_checks: disabled
add_pip_as_python_dependency : false
rollback_enabled: false
notify_outdated_conda: false
```

## How to fix jobs bcbio jobs timing out
The O2 cluster can take a really long time to schedule jobs. If you are having problems with bcbio timing out, set your --timeout parameter to something high, like this:
```bash
/n/app/bcbio/dev/anaconda/bin/bcbio_nextgen.py ../config/bcbio_ensembl.yaml -n 72 -t ipython -s slurm -q short -r --tag feany --timeout 6000 -t 0-11:00
```
