:orphan:

.. _installation:

============
Installation
============

*cctk* requires `Python 3.7+ <python.org>`_, `numpy <numpy.org>`_, and `networkx <networkx.github.io>`_. 
A full list of requirements can be found in ``env.yml`` (`link <https://github.com/ekwan/cctk/blob/master/env.yml>`_.

--------------------------------------------------
Installing with a working Python 3.7+ environment:
--------------------------------------------------

Simply run::

    $ pip install cctk

-----------------------------------------------------
Installing without a working Python 3.7+ environment:
-----------------------------------------------------

If you have a different version of Python (e.g. Python 2.7), you can use a `conda <https://docs.conda.io/en/latest/>`_ environment to run *cctk* without breaking existing packages.

- Install ``conda``/``miniconda``.
- Use ``env.yml`` to create a ``conda`` environment called ``cctk`` and install *cctk*::

    $ conda env create -f env.yml

Now, run ``conda activate cctk`` to enter the *cctk* Python environment (and ``conda deactivate`` to leave). (More complete guides to ``conda`` usage can be found elsewhere.)

---------------
Upgrading cctk
---------------

To get the latest release of cctk, navigate to the correct ``conda`` environment and run::

    $ pip install --upgrade cctk

To install the development version (may be unstable!), run::

    $ pip install --upgrade git+git@github.com:ekwan/cctk.git@master 



