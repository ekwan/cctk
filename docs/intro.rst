.. _intro:

============
Introduction
============

Overview
========

*cctk* is an open-source Python 3 package designed to automate generation and analysis of computational chemistry files. 

Potential uses for *cctk* include:

- Monitoring one or many geometry optimizations.
- Extracting geometry from output files, changing geometric parameters, and creating new input files.
- Calculating molecular properties (e.g. NICS) along a reaction coordinate.
- Screening different functionals and basis sets.
- Generating potential energy surfaces in one or more dimensions (e.g. More O'Ferrall-Jencks plots).

Philosophy
__________

Computational chemistry has become a cornerstone of modern research in all branches of chemistry. 
Most computational chemistry is still performed by manual processing of input/output files, with the aid of 3D structure editors like GaussView, Avogadro, or Molden.
For large projects, the creation and analysis of such files can be prohibitively time-consuming. 

*cctk* provides an object-oriented framework in which files and molecules can be handled simply as Python objects, allowing users to write clean and intuitive code. 
This permits quick generation and analysis of large numbers of jobs, and automates much of the "busy work" assocatied with computational chemistry. 

In general, all *cctk* programs are divided into three parts:

1. Converting computational chemistry files into Python objects.
2. Performing manipulations on these Python objects. 
3. Converting the Python objects back to computational chemistry files, or outputting data in a human-readable format. 

Compatible Filetypes
____________________

Currently *cctk* supports:

- Gaussian 16 ``.out`` (read) and ``.gjf`` (read/write).
- ``.xyz`` (read/write)
- ``.mol2`` (read)
- ``.mae`` (read)
- Orca ``.inp`` (write)

We hope to add more filetypes in future releases. 

Getting Started
==============

Installation
____________

*cctk* requires ``numpy``, ``networkx``, and ``importlib_resources``. A full list of requirements can be found in ``env.yml``.

**Installing with a working Python 3.7+ environment:**

Simply run::

    $ pip install cctk

**Installing without a working Python 3.7+ environment:**

If you have a different version of Python (e.g. Python 2.7), you can use a ``conda`` environment to run *cctk* without breaking existing packages.

1. Install ``conda``/``miniconda``.
2. Use ``env.yml`` to create a Conda environment called ``cctk`` and install *cctk*::

    $ cd cctk
    $ conda env create -f env.yml

Now, run ``conda activate cctk`` to enter the ``cctk`` Python environment (and ``conda deactivate`` to leave). 

(More complete guides to ``conda`` usage can be found elsewhere.)

Using *cctk*
____________

*cctk* comes with several standalone Python scripts which can be used without modification. 

These include: 

- ``analyze.py``, for analysis of multiple Gaussian jobs
- ``monitor.py``, for analysis of a single Gaussian job
- ``resubmit.py``, for resubmitting one or many Gaussian jobs
- ``hammett_swap.py``, for modifying structures through group substitution

These scripts can be found in the ``/scripts`` folder (with appropriate documentation). 

For more sophisticated uses of *cctk*, read on!
