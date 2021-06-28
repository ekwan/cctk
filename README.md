
# cctk

[![PyPI version](https://badge.fury.io/py/cctk.svg)](https://pypi.python.org/pypi/cctk/)
[![Doc status](https://readthedocs.org/projects/pip/badge/)](https://cctk.rtfd.io)
[![Downloads](https://img.shields.io/pypi/dm/cctk.svg)](https://pypi.python.org/pypi/cctk/)
[![Code Quality](https://img.shields.io/lgtm/grade/python/g/ekwan/cctk.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/ekwan/cctk/context:python)

*a Python-based computational chemistry toolkit*

*cctk* simplifies the computational modeling of organic reactions and small molecule structures by automating routine interactions with quantum chemistry software packages:

 - **input file creation**: conformer enumeration, job keyword manipulations, constrained potential energy surface creation
 - **method screening**: creating jobs that screen grids of DFT methods and basis sets
 - **job monitoring**: identification of job status, progress of optimizations, and resubmission of failed jobs
 - **data extraction**: geometries, energies, molecular properties (e.g. charges or NMR shieldings), or geometric parameters (distances, angles, dihedrals) from output files
 - **data analysis**: easy export for statistical analysis or visualization

A quick-start guide is [available](https://cctk.readthedocs.io/en/latest/quick-start.html). More documentation is [here](https://cctk.readthedocs.io/).
 
*cctk* is primarily designed for use with [Gaussian 16](https://gaussian.com). Some support is provided for other file formats (`.xyz`, `.mol2`, `.pdb`, [Schrodinger](https://www.schrodinger.com) `mae`, and [Orca](https://sites.google.com/site/orcainputlibrary/) `.inp`/`.out`).

## Installation

### First Time

*cctk* is easy to install! It should work on any system where Python works.

With Python 3.7 or later, type:

```
pip install cctk
```

If you don't have [pip](https://pypi.org/project/pip/) or virtual environments available on your system, then we recommend installing Anaconda first:

1. Go to [https://www.anaconda.com/distribution/](https://www.anaconda.com/distribution/). Download the Python 3 installer appropriate to your system and run it.

2. Create a virtual environment to use with *cctk*:

 ```
 conda create --name cctk python=3.8
 ```

3. Now activate the virtual environment:

 ```
 conda activate cctk
 ```

To use *cctk*, you will need to place this command at the beginning of your Python scripts:

```
import cctk
```

The [documentation](https://cctk.readthedocs.io/) contains many examples of how to write *cctk* scripts.

### Upgrading

*cctk* is undergoing active development. To upgrade to the latest stable release:

```
pip install --upgrade cctk
```

To install the development version, which may be unstable, run:

```
$ pip install --upgrade git+git@github.com:ekwan/cctk.git@master 
```

Alternatively, clone the repository. Then, from within the repository folder, run:

```
pip install --upgrade .
```

### Building Documentation

If you want to read the *cctk* documentation locally, you can build it by going to the `docs` folder and typing:

```
make html
```

This command will require the `sphinx` and `sphinx-bootstrap-theme` packages to be installed first. Once generated, the documentation will be available locally at: `docs/_build/html/index.html`.

## Fine Print

### Package Details 

- `cctk/` contains the Python modules for *cctk* and the accompanying static data files.  
- `docs/` contains the code needed to generate the documentation.  
- `scripts/` contains pre-defined scripts that use *cctk* to quickly analyze and manipulate one or many output files.  
- `test/` contains code to test *cctk* and accompanying files.  
- `tutorial/` contains detailed tutorials on how to use *cctk* on complex, real-world problems.  

*cctk* requires Python 3.7+, [`numpy`](https://numpy.org/), and [`networkx`](https://networkx.github.io/).
A full list of requirements can be found in `env.yml`. 

### External Data:

*cctk* depends on some external data (`cctk/data/`):

- Atomic weights are taken from the 
[NIST website](https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&all=all&ascii=ascii2&isotype=some) 
and stored in `cctk/data/isotopes.csv`.
- Covalent radii are taken from 
[*Dalton Trans.* **2008**, 2832](https://pubs.rsc.org/en/content/articlelanding/2008/dt/b801115j#!divAbstract) 
and stored in `cctk/data/covalent_radii.csv`.
(When multiple atomic types were specified, the one with longer bond distances was adopted for simplicity).
- van der Waals radii were taken from
[*J. Am. Chem. Soc.* **1983**, *105*, 5220](https://pubs.acs.org/doi/10.1021/ja00354a007), 
[*Inorg. Mater.* **2001**, *37*, 871](https://link.springer.com/article/10.1023/A:1011625728803), and
[*J. Phys. Chem. A*, **2009**, *113*, 5806](https://pubs.acs.org/doi/10.1021/jp8111556) and stored in `cctk/data/vdw_radii.csv`.

### Authors and Community:

*cctk* is an ongoing project led by Corin Wagen and Eugene Kwan.

There is a Slack workspace for *cctk* users and developers: please email ``cwagen@g.harvard.edu`` or ``ekwan16@gmail.com`` for access.

### How to Cite:

Wagen, C.C.; Kwan, E.E. *cctk* **2020**, [www.github.com/ekwan/cctk](https://www.github.com/ekwan/cctk).

### License:

This project is licensed under the Apache License, Version 2.0.  Please see `LICENSE` for full terms and conditions. 

*Copyright 2020 by Corin Wagen and Eugene Kwan*
