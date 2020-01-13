# cctk
## Computational Chemistry Toolkit

*This is a Python 3-based library for working with computational chemistry data*.

## Contents: 
 - [Overview](#overview) 
 - [Requirements](#requirements)
 - [Installation](#installation)
 - [Structure](#structure)
 - [Documentation](#requirements)
 - [Technical Details](#technical-details)

## Overview:

*cctk* is an open-source Python package designed to automate generation and analysis of computational chemistry files. 

Potential uses for *cctk* include: 
 - Monitoring one or many geometry optimizations. 
 - Extracting geometry from output files, changing geometric parameters, and creating new input files. 
 - Calculating molecular properties (e.g. NICS) along a reaction coordinate. 
 - Screening different functionals and basis sets. 
 - Generating potential energy surfaces in one or more dimensions (e.g. More O'Ferrall-Jencks plots). 
 
 For hands-on examples of how *cctk* can be used, 
 refer to the [tutorials](https://github.com/ekwan/cctk/tree/master/tutorial). 
 
### Compatible File Types:
 - Gaussian 16 `.out` (read) and `.gjf` (read/write).
 - `.xyz` (read/write)
 - `.mol2` (read)
 - `.mae` (read)

## Installation:

*cctk* requires Python 3.7, [`numpy`](https://numpy.org/), and [`networkx`](https://networkx.github.io/).
A full list of requirements can be found in `requirements.txt`.

The preferred setup method is as follows: 

1. Install Python 3.7+ or `conda`/`miniconda`


## Structure: 

Most *cctk* programs follow a rough 3-part outline: 

1. Read in data from an output file (or files). 
1. Perform some transformation (adding an atom, changing bond lengths, etc.). 
1. Output an input file. 

## Requirements:
* Python 3.7 or later
* Numpy
* NetworkX
* Sphinx

Use of a package management system like `conda` or `miniconda` is recommended. To install all requirements, run:

```
pip install -r requirements.txt
```

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
