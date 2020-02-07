# cctk
## Computational Chemistry Toolkit

*This is a Python 3-based library for working with computational chemistry data*.

## Contents: 
 - [Overview](#overview) 
 - [Installation](#installation)
 - [Contents](#contents)
 - [Documentation](#documentation)
 - [Technical Details](#technical-details)
 - [Authors](#authors)
 - [How to Cite](#how-to-cite)
 - [License](#license)

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

*cctk* requires Python 3.7+, [`numpy`](https://numpy.org/), and [`networkx`](https://networkx.github.io/).
A full list of requirements can be found in `env.yml`.

#### Installing with a working Python 3.7+ environment:

Simply run: 
```
$ pip install cctk
```

#### Installing without a working Python 3.7+ environment:

If you have a different version of Python (e.g. Python 2.7), 
you can use a `conda` environment to run *cctk* without breaking existing packages.

1. Install [`conda`](https://docs.conda.io/en/latest/)/[`miniconda`](https://docs.conda.io/en/latest/miniconda.html).

2. Use `env.yml` to create a Conda environment called `cctk` and install *cctk*:

```
$ cd cctk
$ conda env create -f env.yml
```

Now, run `conda activate cctk` to enter the *cctk* Python environment (and `conda deactivate` to leave).
(More complete guides to `conda` usage can be found elsewhere.)

## Contents: 

- `cctk/` contains the Python modules for *cctk* and the accompanying static data files.  
- `docs/` contains the code needed to generate the documentation.  
- `scripts/` contains pre-defined scripts that use *cctk* to quickly analyze and manipulate one or many output files.  
- `test/` contains code to test *cctk* and accompanying files.  
- `tutorial/` contains detailed tutorials on how to use *cctk* on complex, real-world problems.  

## Documentation:

The documentation for *cctk* can be found on [Read the Docs](https://cctk.readthedocs.io). 

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

## Authors:

*cctk* was written by Corin Wagen and Eugene Kwan at the Department of Chemistry and Chemical Biology at Harvard University. 
Please email `cwagen@g.harvard.edu` with any questions or bug reports; we will do our best!

## How to Cite:

Wagen, C.C.; Kwan, E.E. *cctk* **2020**, [www.github.com/ekwan/cctk](https://www.github.com/ekwan/cctk).

## License:

This project is licensed under the Apache License, Version 2.0: see `LICENSE` for full terms and conditions. 

*Copyright 2020 by Corin Wagen and Eugene Kwan*
