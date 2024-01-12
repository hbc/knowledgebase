## Conda

Every system has a Python installation, but you don't necessarily want to use that. Why not? That version is typically outdated and configured to support system functions. Most tools require specific versions of Python and depedencies so you need more flexibility.

**Solution?**

Set up a full-stack scientific Python deployment is to **use a Python distribution** (Anaconda or Miniconda). It is an installation of Python with a set of curated packages which are guaranteed to work together.
 

## Setting up Python distribution on O2

You can install it to your home directory, though not needed as O2 has a miniconda module available. 

By default, miniconda and conda envs are installed under user home space. 

### Conda Environments
Environments allow you to creata an isolated, reproducible environments where you have fine-tuned control over python version, all packages and configuration. _This is always recommended over using the default environment._

To create an environment using Pytho 3.9 and the numpy package:

```bash
$ conda create --name my_environment python=3.9 numpy
```

Now that you have created it, you need to activate it. Once activated, all installations of tools are specific to that environment; do this using `conda install`. It is a configured space where you can run analyses reproducibly.

```bash
$ conda activate my_environment
```

When you are done you can deactivate the environment or close it:

```bash
conda deactivate
```

The environments and associated libs are located in `~/home/minconda3/envs`. Both miniconda and the created environments can occupy a lot of space and max out your home directory! 

**Solution: Create the conda env in another space**

For conda envs, you can use the full path outside of home when creating env:

```bash
module purge
module load miniconda3/23.1.0
conda create -p /path/to/somewhere/not/home/myEnv -c channel1 packagex==1.23
```

> **NOTE:** It's common that installing packages using Conda is slow or fails because Conda is unable to resolve dependencies. To get around this, we suggest the use of Mamba.

**Installing lots of dependency packages?**

You can do this easily by creating a YAML file, for example 

### Channels 

> Adapted from [An Introduction to Earth and Environmental Data Science](https://earth-env-data-science.github.io/lectures/environment/python_environments.html) 
