## Conda

Every system has a Python installation, but you don't necessarily want to use that. Why not? That version is typically outdated and configured to support system functions. Most tools require specific versions of Python and depedencies so you need more flexibility.

**Solution?**

Set up a full-stack scientific Python deployment **using a Python distribution** (Anaconda or Miniconda). It is an installation of Python with a set of curated packages which are guaranteed to work together.
 

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
conda create -p /path/to/somewhere/not/home/myEnv python=3.9 numpy
```

> **NOTE:** It's common that installing packages using Conda is slow or fails because Conda is unable to resolve dependencies. To get around this, we suggest the use of Mamba.

**Installing lots of dependency packages?**

You can do this easily by creating a yaml file, for example `environment.yaml` below was used to install Pytables:

```bash
name: pytables
channels:
  - defaults
dependencies:
  - python=3.9*
  - numpy >= 1.19.0
  - zlib
  - cython >= 0.29.32
  - hdf5=1.14.0
  - numexpr >= 2.6.2
  - packaging
  - py-cpuinfo
  - python-blosc2 >= 2.3.0
```

Now to create the environment we reference the file in the command:

```bash
conda env create -f environment.yaml
```

### Channels 

Where do conda packages come from? The packages are hosted on conda “channels”. From the conda pages:

_"Conda channels are the locations where packages are stored. They serve as the base for hosting and managing packages. Conda packages are downloaded from remote channels, which are URLs to directories containing conda packages. The conda command searches a set of channels."_

Using `-c` you can specify whihc channels you want conda to search in for packages.

> Adapted from [An Introduction to Earth and Environmental Data Science](https://earth-env-data-science.github.io/lectures/environment/python_environments.html) 
