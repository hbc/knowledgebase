## How to fix jobs bcbio jobs timing out

The O2 cluster can take a really long time to schedule jobs. If you are having problems with bcbio timing out, set your --timeout parameter to something high, like this:

eg.  /n/app/bcbio/dev/anaconda/bin/bcbio_nextgen.py ../config/bcbio_ensembl.yaml -n 72 -t ipython -s slurm -q short -r --tag feany --timeout 6000 -t 0-11:00

