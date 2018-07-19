The O2 cluster can take a really long time to schedule jobs. If you are having problems with bcbio timing out, set your time limit (-t) to something very high like 6000 minutes
eg.  /n/app/bcbio/dev/anaconda/bin/bcbio_nextgen.py ../config/bcbio_ensembl.yaml -n 72 -t ipython -s slurm -q short -r **t=0-6000:00** --tag feany
