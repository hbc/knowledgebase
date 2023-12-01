## Setting up conda environments on O2
Conda is a package management system which can be used to create Python and R development environments on many different platforms. 

You can install it to your home directory, though not needed as O2 has a minoconda module. 

By default, miniconda and conda envs are installed under user home space. 

Both miniconda and the created environments can occupy a lot of space and max out your home directory! 

**Solution: Create the conda env in another space**

For conda envs, you can use full path outside of home when creating env:

```bash
module purge
module load miniconda3/23.1.0
mamba create -p /path/to/somewhere/not/home/myEnv -c channel1 packagex==1.23
```

Example for package install for [dictys](https://github.com/pinellolab/dictys/tree/master):

```
conda create -p  /n/data1/cores/bcbio/meeta_mistry/conda/dictys -y -c conda-forge python=3.9 mamba
``
 
