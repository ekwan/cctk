.. _intro:

====
cctk
====

computational chemistry toolkit
====

Overview:
____

*cctk* is an open-source Python 3 package designed to automate generation and analysis of computational chemistry files. 

Potential uses for *cctk* include:

- Monitoring one or many geometry optimizations.
- Extracting geometry from output files, changing geometric parameters, and creating new input files.
- Calculating molecular properties (e.g. NICS) along a reaction coordinate.
- Screening different functionals and basis sets.
- Generating potential energy surfaces in one or more dimensions (e.g. More O'Ferrall-Jencks plots).

Compatible Filetypes:
____

Currently *cctk* supports:

- Gaussian 16 ``.out`` (read) and ``.gjf`` (read/write).
- ``.xyz`` (read/write)
- ``.mol2`` (read)
- ``.mae`` (read)
- Orca ``.inp`` (write)

We hope to add more filetypes in future releases. 

Installation:
____

*cctk* requires ``numpy``, ``networkx``, and ``importlib_resources``. A full list of requirements can be found in ``env.yml``.

**Installing with a working Python 3.7+ environment:**

Simply run:

``$ pip install cctk``

**Installing without a working Python 3.7+ environment:**

If you have a different version of Python (e.g. Python 2.7), you can use a ``conda`` environment to run *cctk* without breaking existing packages.

1. Install ``conda``/``miniconda``.
2. Use ``env.yml`` to create a Conda environment called ``cctk`` and install *cctk*:

``$ cd cctk``
``$ conda env create -f env.yml``

Now, run ``conda activate cctk`` to enter the ``cctk`` Python environment (and ``conda deactivate`` to leave). (More complete guides to ``conda`` usage can be found elsewhere.)

Authors:
____
*cctk* was written by Corin Wagen and Eugene Kwan at the Department of Chemistry and Chemical Biology at Harvard University. 
Please email ``cwagen@g.harvard.edu`` with any questions or bug reports; we will do our best!

How to Cite:
____
Wagen, C.C.; Kwan, E.E. *cctk* 2020, www.github.com/ekwan/cctk.

License:
____
This project is licensed under the Apache License, Version 2.0.

*Copyright 2020 by Corin Wagen and Eugene Kwan*
