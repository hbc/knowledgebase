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
mkdir -p ${HOME}/local/share/galaxy
ln -s /n/app/bcbio/biodata/galaxy/tool-data ${HOME}/local/share/galaxy/tool-data
export PATH="${HOME}/local/bin:$PATH"
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

## How to run a one-node bcbio job (multicore, not multinode)
it just runs a bcbio job on one node of the cluster (no IPython)
[More slurm options](https://wiki.rc.hms.harvard.edu/display/O2/Using+Slurm+Basic#UsingSlurmBasic-sbatchoptionsquickreference)

```
#!/bin/bash

# https://slurm.schedmd.com/sbatch.html

#SBATCH --partition=priority        # Partition (queue)
#SBATCH --time=3-00:00              # Runtime in D-HH:MM format
#SBATCH --job-name=bcbio            # Job name - any name
#SBATCH -c 10                       # cores per task 
#SBATCH --mem-per-cpu=10G           # Memory needed per CPU or use --mem to limit total memory
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL             # Type of email notification (BEGIN, END, FAIL, ALL) by default goes to the email associated with O2 accounts
#SBATCH --mail-user=abc123@hms.harvard.edu   # Email to which notifications will be sent

bcbio_nextgen.py ../config/illumina_rnaseq.yaml -n 10
```
