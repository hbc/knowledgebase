# Jupyter notebooks

This post is for those who are interested in running notebooks seemingly on O2. There is a well written documentation about running jupyter notebook on O2. you can find it here. https://wiki.rc.hms.harvard.edu/display/O2/Jupyter+on+O2. However, this involves multiple steps, opening bunch of terminals at times, and importantly finding an unused port every-time. I found it quite cumbersome and annoying, so i spent some time solving it. It took me a while to nail it down, with the help of FAC RC, but they suggested a simpler solution. If you wish to run jupyter/R notebook on O2 {where your data sits},
Here what you need to do : 

Install https://github.com/aaronkollasch/jupyter-o2 by `pip install jupyter-o2` on your local terminal. 

Run `jupyter-o2 --generate-config` on command line. 
This will generate the configuration file and will tell you where it is located. Un comment the fields which are not needed. Since configuration file is the key, find attach the a template for use, you need to change your credentials though.

You are all set to run notebook on O2 from your local machine now, without logging into server.  
Now at your local terminal Run `jupyter-o2 notebook` for python notebooks. Alternatively you can also do `jupyter-o2 lab` for R/python
This will ask you a paraphrase, you should enter your ecommons password as paraphrase. 
Boom!!! you are good to go! Happy Pythoning :):)  
If you wish you run R notebooks on O2, refer this. https://docs.anaconda.com/anaconda/navigator/tutorials/r-lang/


# Example code

Just to add, in the HMS-RC documentation they suggested any ports over 50000. To give examples of logging into a jupyter notebook session I have provided the code below.

## Creating a Jupyter notebook

Log onto a login node

```
# Log onto O2 using a specific port - I used '50000' in this instance - you can choose a different  port and just replace the 50000 with the number of your specific port
ssh -Y -L 50000:127.0.0.1:50000 ecommons_id@o2.hms.harvard.edu 
```

Once on the login node, you can start an interactive session specifying the port with `--tunnel`

```
# Create interactive session
srun --pty -p interactive -t 0-12:00 --x11 --mem 128G --tunnel 50000:50000 /bin/bash
```

Load the modules that you will need

```
# Load modules
module load gcc/9.2.0 python/3.8.12
```

Create environment for running analysis (example here is for velocity)

```
# Create virtual environment (only do this once)
virtualenv velocyto --system-site-packages
```

Activate virtual environment

```
# Activate virtual environment
source velocyto/bin/activate
```

Install Jupyter notebook and any other libraries (only need to do this once)

```
# Install juypter notebook
pip3 install jupyter

# Install any other libraries needed for analysis (this is for velocity)
pip3 install numpy scipy cython numba matplotlib scikit-learn h5py click
pip3 install velocyto
pip3 install scvelo
```

To create a Jupyter notebook run the following (again instead of 50000, use your port #):

```
# Start jupyter notebook
jupyter notebook --port=50000 --browser='none'
```

## Logging onto an existing notebook

```
# Log onto O2 using a specific port - I used '50000' in this instance - you can choose a different  port and just replace the 50000 with the number of your specific port
ssh -Y -L 50000:127.0.0.1:50000 ecommons_id@o2.hms.harvard.edu 

# Create interactive session
srun --pty -p interactive -t 0-12:00 --x11 --mem 128G --tunnel 50000:50000 /bin/bash

# Load modules
module load gcc/9.2.0 python/3.8.12

# Activate virtual environment
source velocyto/bin/activate

# Open existing notebook
jupyter notebook name_of_notebook.ipynb --port=PORT --browser='none'
