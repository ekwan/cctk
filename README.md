# cctk
## Computational Chemistry Toolkit

*This is a Python 3-based library for working with computational chemistry data*.

## Contents: 
 - [Overview](#overview) 
 - [Installation](#installation)
 - [Contents](#contents)
 - [Documentation](#documentation)
 - [Technical Details](#technical-details)

## Overview:

*cctk* is an open-source Python package designed to automate generation and analysis of computational chemistry files. 

Potential uses for *cctk* include: 
 - Monitoring one or many geometry optimizations. 
 - Extracting geometry from output files, changing geometric parameters, and creating new input files. 
 - Calculating molecular properties (e.g. NICS) along a reaction coordinate. 
 - Screening different functionals and basis sets. 
 - Generating potential energy surfaces in one or more dimensions (e.g. More O'Ferrall-Jencks plots). 
 
 For examples of how *cctk* can be used, 
 refer to the [tutorials](https://github.com/ekwan/cctk/tree/master/tutorial). 
 
### Compatible File Types:
 - Gaussian 16 `.out` (read) and `.gjf` (read/write).
 - `.xyz` (read/write)
 - `.mol2` (read)
 - `.mae` (read)
 - Orca `.inp` (write)

## Installation:

*cctk* requires Python 3.7, [`numpy`](https://numpy.org/), and [`networkx`](https://networkx.github.io/).
A full list of requirements can be found in `requirements.txt`.

The preferred setup method is as follows: 

1. Install Python 3.7+ and 
[`conda`](https://docs.conda.io/en/latest/)/[`miniconda`](https://docs.conda.io/en/latest/miniconda.html)
2. Create a new environment called `cctk` (the name is arbitrary): 

```
$ conda create --name cctk python=3.8
$ source activate cctk
```

3. `git clone` this repository, or download the `.zip` file, 
and use `pip` to install the required packages into the `cctk` environment:

```
$ git clone git@github.com:ekwan/cctk.git
$ cd cctk
$ pip install -r requirements.txt
```

4. Add *cctk* to the `PYTHONPATH` in your bash configuration file (`~/.bashrc`) by adding the following line:

```
export PYTHONPATH="$PYTHONPATH:/path/to/cctk/"
```
(be sure to replace `/path/to/cctk/` with whatever's correct for your system!)

5. Restart bash (or type `$ source ~/.bashrc)` to allow these changes to take effect. 

You should now be able to import *cctk* as a Python library anywhere on your system. 


## Contents: 

- `cctk/` contains the Python modules for *cctk* and the accompanying static data files.  
- `docs/` contains the code needed to generate the documentation.  
- `scripts/` contains pre-defined scripts that use *cctk* to quickly analyze and manipulate one or many output files.  
- `test/` contains code to test *cctk* and accompanying files.  
- `tutorial/` contains detailed tutorials on how to use *cctk* on complex, real-world problems.  

## Documentation:

To build the documentation, run: 

```
cd docs/
sphinx-apidoc -o . ../cctk/
make html
```

The documentation files can then be found in `docs/_build/html`.

## Technical Details: 

### External Data:

*cctk* depends on some external data, stored in `cctk/data/`:
- Atomic weights are taken from the 
[NIST website](https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&all=all&ascii=ascii2&isotype=some) 
and stored in `cctk/data/isotopes.csv`.
- Covalent radii are taken from 
[**Dalton Trans.** *2008*, 2832&ndash;2838](https://pubs.rsc.org/en/content/articlelanding/2008/dt/b801115j#!divAbstract) 
and stored in `cctk/data/covalent_radii.csv`.
(When multiple atomic types were specified, the one with longer bond distances was adopted for simplicity).

*Written by Eugene Kwan and Corin Wagen.*
