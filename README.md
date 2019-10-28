# cctk
## Computational Chemistry Toolkit

*This is a Python 3-based library for working with computational chemistry data*.

#### Requirements:
* Python 3.x
* Numpy
* NetworkX
* Sphinx

Use of a package management system like `conda` or `miniconda` is heavily recommended. To install all requirements, run:

```
pip install -r requirements.txt
```

#### Documentation

To build the documentation, run: 

```
cd docs/
make html
```

#### External Data:

Atomic weights taken from the [NIST website](https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&all=all&ascii=ascii2&isotype=some). 

Covalent radii taken from [**Dalton Trans.** *2008*, 2832&ndash;2838](https://pubs.rsc.org/en/content/articlelanding/2008/dt/b801115j#!divAbstract). (when multiple atomic types were specified, the one with longer bond distances was adopted).

*Written by Eugene Kwan and Corin Wagen.*
