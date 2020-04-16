:orphan:

.. _overview: 
.. |br| raw:: html

  <br/>

.. |nbsp| raw:: html 

   &nbsp;

========
Overview
========

`cctk <https://www.github.com/ekwan/cctk>`_: a Python-based computational chemistry toolkit.

*cctk* simplifies routine tasks in computational chemistry: preparing input files with scripts,
checking whether jobs ran successfully, extracting energies and geometries, etc. All *cctk*
operations are carried out using Python scripts.  The prototypical workflow involves:

1. Reading in output files from a quantum chemistry program like Gaussian.
2. Analyzing the extracted data (e.g., determining which structure is lowest
   in energy).
3. Writing out new input files for further calculations.
4. Further analysis or visualization with
   `pandas <https://https://pandas.pydata.org/>`_ or
   `matplotlib <https://matplotlib.org/>`_.

==============
*cctk* Objects
==============

Use these three main classes to interact with external quantum chemistry programs:

""""""""""""""""
1.  ``Molecule``
""""""""""""""""
    A single molecular geometry.

.. rst-class:: full-width-table

    =================================   ===========================================
    Field                               Description
    =================================   ===========================================
    ``molecule.atomic_numbers``         the atomic number for each atom
    ``molecule.geometry``               xyz coordinates
    ``molecule.bonds``                  the connectivity as a `networkx <https://networkx.github.io>`_ graph
    ``molecule.charge``                 the overall charge
    ``molecule.multiplicity``           the spin multiplicity
    =================================   ===========================================

|nbsp|
    All arrays that refer to atoms in *cctk* are 1-indexed (i.e., 1, 2, ..., n).
    Thus, both the ``atomic_numbers`` and ``geometry`` fields are 1-indexed.
    In contrast, all arrays that refer to non-atoms are 0-indexed.

    Various methods are available to measure or set geometric parameters (bond distances, bond angles,
    or dihedral angles).

""""""""""""""""
2.  ``Ensemble``
""""""""""""""""
    A collection of molecules and associated properties.

    Each ``Molecule`` in the ``Ensemble`` is associated with its properties
    (filenames, energies, NMR shieldings, etc.) using a dictionary.  For example,
    a conformation of pentane might be mapped to this ``dict``::

        properties_dict = {
             'energy': -0.0552410743198,
             'scf_iterations': 2,
             'link1_idx': 0,
             'filename': 'test/static/pentane_conformation_1.out',
             ... }
    
 
    To access ``Ensemble`` information, use the following syntax:
    
.. rst-class:: full-width-table
    
    ==========================================      ==============================================================================
    Syntax                                          Result
    ==========================================      ==============================================================================
    ``ensemble.molecules``                          iterator over all molecules
    ``ensemble.molecules[i]``                       the *i*-th molecule (0-indexed)
    ``ensemble.molecules[1:3]``                     the second and third molecules as a list
    ``ensemble.molecules[-1]``                      the last molecule
    ``ensemble.items()``                            iterator over (molecule, property dictionary) tuples
    ``ensemble.get_properties_dict(molecule)``      the property dictionary associated with ``molecule`` 
    ``ensemble[:,"energy"]``                        one-dimensional array of energies, with ``None`` as a placeholder for any missing data
    ``ensemble[:,["filename","energy"]]``           two-dimensional array of filenames and energies, with ``None`` as a placeholder for any missing data
    ``ensemble.molecule_list()``                    list of molecules
    ``ensemble.properties_list()``                  list of the property dictionaries
    ``ensemble[0]``                                 ``Ensemble`` containing the first molecule and its properties
    ``ensemble[0:2]``                               ``Ensemble`` containing the first and second molecules and their properties
    ==========================================      ==============================================================================

|nbsp|
    Thus, Ensembles can be indexed or sliced to return smaller Ensembles.  Note that while all
    such sub-Ensembles are new ``Ensemble`` objectes, they are essentially views of the original
    ``Ensemble``, rather than deep copies.

    A ``ConformationalEnsemble`` is a special case of an ``Ensemble`` in which each structure
    corresponds to the same molecule.  This allows for RMSD calculation, structural alignment,
    and redundant conformer elimination to be carried out as desired (see tutorials).

"""""""""""""""""""
3. ``GaussianFile``
"""""""""""""""""""
    The results of a Gaussian job or the contents of an input file::
    
		gaussian_file = cctk.GaussianFile.read_file(filename)
    
    ``filename`` may be a Gaussian output file (``.out``/``.log``) or a Gaussian input file
    (``.gjf``/``.com``).

    **Important:** *cctk* assumes that all Gaussian jobs will be run in verbose mode (``#p``
    in the route card).  **Parsing will not work correctly without** ``#p`` **.**

    As usual, molecules and their properties are stored in ``gaussian_file.ensemble``::

        ensemble = first_link.ensemble
        energies = list(ensemble[:,"energy"])
        # [-40.5169484082, -40.5183831835, -40.5183831835])
        
        ensemble = second_link.ensemble
        shieldings = ensemble[-1,"isotropic_shielding"]
        # [192.9242, 31.8851, 31.8851, 31.8851, 31.8851]
    
    Per *cctk* convention (*vide infra*), ``energies`` is 0-indexed, but ``shieldings`` is
    1-indexed.  (The ``-1`` refers to the last geometry.)
    
    (*Note:* if a Gaussian input file is read, no properties will be available, so
    the `properties_dict` for each molecule will be empty.)

    Some Gaussian output files are composites of multiple jobs using the
    `Link1 <http://gaussian.com/input/>`_ directive.  In that case,
    ``GaussianFile.read_file(filename)`` will return one ``GaussianFile``
    object per Link1 section.

    For example, this is a two-step job::

        gaussian_file = cctk.GaussianFile.read_file("test/static/methane2.out")
        assert len(gaussian_file), 2
        first_link = gaussian_file[0]
        second_link = gaussian_file[1]

    *cctk* will also interpret common job types via the ``cctk.JobType`` enum::

        # first_link.job_types = [JobType.OPT, JobType.FREQ, JobType.SP]

.. rst-class:: full-width-table

    =================================   ===========================================
    Field                               Description
    =================================   ===========================================
    ``gaussian_file.ensemble``          ``Ensemble`` containing intermediate geometries and molecular properties
    ``gaussian_file.job_types``         list of what kind of jobs were run
    ``gaussian_file.success``           number of successful terminations
    ``gaussian_file.link0``             dictionary containing `Link0 <https://gaussian.com/link0/>`_ information (memory, processors, checkpoint filename, etc.)
    ``gaussian_file.route_card``        route card (must start with ``#p``)
    ``gaussian_file.title``             title of Gaussian file
    ``gaussian_file.footer``            footer (optional)
    =================================   ===========================================

    |nbsp|

    Limited support for other file formats is available (see Features section of documentation).

========
Indexing
========

In *cctk*, **arrays whose contents refer to atoms are always 1-indexed; other arrays are 0-indexed.**

Thus, arrays of atomic numbers, positions, or NMR shieldings are 1-indexed, while arrays of
molecules, files, or molecular property values are 0-indexed.

1-indexed arrays are implemented via ``cctk.OneIndexedArray``, a custom subclass of ``np.ndarray``.  
For example::

    molecule.geometry[1]

will return the coordinates of the first atom of the ``Molecule``.  However::

    ensemble.molecules[0]

returns the first molecule of the ``Ensemble``.


