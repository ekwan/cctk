:orphan:

.. _installation:

============
Installation
============

----------
First Time
----------

*cctk* is easy to install! It should work on any system where Python works.

With Python 3.7 or later, type::

    pip install cctk

If you don't have ``pip`` or virtual environments available on your system, then we recommend installing Anaconda first:

1. Go to https://www.anaconda.com/distribution/. Download the Python 3 installer appropriate to your system and run it.

2. Create a virtual environment to use with *cctk*::

    conda create --name cctk python=3.8

3. Now activate the virtual environment::

    conda activate cctk

To use *cctk*, you will need to place this command at the beginning of your Python scripts::

    import cctk

The documentation contains many examples of how to write *cctk* scripts.

---------
Upgrading
---------

*cctk* is undergoing active development. To upgrade to the latest stable release::

    pip install --upgrade cctk

To install the development version, which may be unstable, run::

    pip install --upgrade git+git@github.com:ekwan/cctk.git@master

Alternatively, clone the repository. Then, from within the repository folder, run::

    pip install --upgrade .

----------------------
Building Documentation
----------------------

If you want to read the *cctk* documentation locally, you can build it by going to the docs folder and typing::

    make html

This command will require the ``sphinx`` and ``sphinx-bootstrap-theme`` packages to be installed first. 
Once generated, the documentation will be available locally at: ``docs/_build/html/index.html.``

