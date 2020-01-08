# cctk
## Computational Chemistry Toolkit

*This is a Python 3-based library for working with computational chemistry data*.

## Contents: 
 - [Overview](#overview) 
 - [Requirements](#requirements)
 - [Installation](#installation)
 - [Structure](#structure)
 - [Documentation](#requirements)
 - [External Data](#external-data)

## Overview:

*cctk* is an open-source Python package designed to automate routine computational chemistry tasks. 

Potential uses for *cctk* include: 
 - Monitoring one or many geometry optimizations. 
 - Extracting geometry from output files, changing geometric parameters, and creating new input files. 
 - Calculating molecular properties (e.g. NICS) along a reaction coordinate. 
 - Screening different functionals and basis sets. 
 - Generating potential energy surfaces in one or more dimensions (e.g. More O'Ferrall-Jencks plots). 
 
### Compatible File Types:
 - Gaussian 16 `.out` (read) and `.gjf` (read/write).
 - `.xyz` (read/write)
 - `.mol2` (read)
 - `.mae` (read)

## Installation:

For now, just use `git clone <git url>` - a more advanced way is coming.

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

## External Data:

Atomic weights taken from the [NIST website](https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&all=all&ascii=ascii2&isotype=some). 

Covalent radii taken from [**Dalton Trans.** *2008*, 2832&ndash;2838](https://pubs.rsc.org/en/content/articlelanding/2008/dt/b801115j#!divAbstract). (when multiple atomic types were specified, the one with longer bond distances was adopted).

*Written by Eugene Kwan and Corin Wagen.*
