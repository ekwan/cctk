.. _overview: 

========
Overview
========

`cctk <https://www.github.com/ekwan/cctk>`_ is a Python-based computational chemistry toolkit.
*cctk* simplifies routine tasks in computational chemistry: preparing input files with scripts,
checking whether jobs completed correctly, extracting energies and geometries, etc.  Here, we
will discuss the basic philosophy behind *cctk*.

All *cctk* operations are carried out using Python scripts.  The prototypical workflow
involves:

1. Reading in output files from a quantum chemistry program like Gaussian.
2. Analyzing the extracted data (e.g., determining which structure is lowest
   in energy).
3. Writing out new input files for further calculations.
4. We could also export data to other Python packages like
   `pandas <https://https://pandas.pydata.org/>`_ or
   `matplotlib <https://matplotlib.org/>`_ for further analysis or visualization.

--------------
*cctk* Objects
--------------

Use these four main classes to interact with external quantum chemistry programs:

""""""""""""""""
1.  ``Molecule``
""""""""""""""""
    A single molecular geometry.

    ``molecule.atomic_numbers``: the element for each atom (1-indexed, see below)
    ``molecule.geometry``: the xyz coordinates (also 1-indexed)

""""""""""""""""
2.  ``Ensemble``
""""""""""""""""
    A collection of Molecules. To access:

    =================================   ===========================================
    Syntax                              Result
    =================================   ===========================================
    ``ensemble.molecules[i]``           the *i*-th molecule (0-indexed)
    ``ensemble.molecules[1:3]``         the second and third molecules as a list
    ``ensemble.molecules[-1]``          the last molecule
    ``for m in ensemble.molecules:``    to iterate
    =================================   ===========================================

    Each ``Molecule`` in the ``Ensemble`` is associated with a dictionary of properties (filenames, energies, NMR shieldings, etc.):

    ``ensemble.get_property_dict(molecule)``

    Alternatively, we can iterate:

    ``for molecule,property_dict in ensemble.items():``
    
    A convenient tuple syntax to obtain properties:

    ``ensemble[:,"filename"]``: all filenames as a list
    ``ensemble[:,["filename","energy"]]``: two-dimensional array, with ``None`` as a placeholder for any missing data

    Ensembles can be indexed or sliced to return smaller Ensembles:

    ``ensemble[0]``: an ``Ensemble`` containing only the first molecule and its properties
    ``ensemble[0:2]``: an ``Ensemble`` containing the first and second molecules and their properties

    Note that while all such sub-Ensembles are new ``Ensemble`` objectes, they are essentially views of the original ``Ensemble``, rather than deep copies.

"""""""""""""""""""""""""""""
3. ``ConformationalEnsemble``
"""""""""""""""""""""""""""""
    An ``Ensemble`` where each structure corresponds to the same molecule.  Methods are provided for structural alignment and redundant conformer elimination (see tutorials).

"""""""""""""""""""
4. ``GaussianFile``
"""""""""""""""""""
    The results of a Gaussian job or the contents of an input file.  To create:

    ``cctk.GaussianFile.read_file(filename)``

    To access the contents:

    ``gaussian_file.ensemble``

    Note that this class can also read and write Gaussian input files.  Some limited support for other file formats is available (see tutorials).

--------
Indexing
--------

The general rule in *cctk* is that atoms are 1-indexed and everything else is 0-indexed.
1-indexed arrays are implemented via ``cctk.OneIndexedArray``.  These work just like
``np.ndarray``.  For example::

    molecule.geometry[1]

will return the coordinates of the first atom.  However::

    ensemble.molecules[0]

is needed to get the first molecule.

-------
Example
-------

Here is a simple script for reading in some conformations of pentane.  First, we will
make a few imports::

    import numpy as np
    import cctk
    import glob as glob
    import pandas as pd
    from pandas import DataFrame

Next, we'll figure out which files we want to open::

    path = "test/static/pentane_conformation*.out"
    filenames = sorted(glob.glob(path))
    
Let's read these files into a ``ConformationalEnsemble``::

    conformational_ensemble = cctk.ConformationalEnsemble()
    for filename in filenames:
        gaussian_file = cctk.GaussianFile.read_file(filename)
        ensemble = gaussian_file.ensemble
        molecule = ensemble.molecules[-1]
        property_dict = ensemble.get_property_dict(molecule)
        conformational_ensemble.add_molecule(molecule,property_dict)

These are geometry optimization jobs, so each ``GaussianFile`` contains
an ``Ensemble`` with each geometry step.  Calling `ensemble.molecules[-1]`
provides the last geometry.  Each ``Molecule`` is associated with a property
dictionary::

    {
     'energy': -0.0552410743198,
     'scf_iterations': 2,
     'link1_idx': 0,
     'filename': 'test/static/pentane_conformation_1.out',
     'rms_force': 4.4e-05,
     'rms_displacement': 0.000319,
     'enthalpy': 0.106416,
     'gibbs_free_energy': 0.068028,
     'frequencies': [101.5041, 117.3291, 192.5335, 201.8222, 231.7895, 463.1763, 465.449, 717.7345, 778.6405, 876.373, 915.2653, 972.8192, 974.4666, 1071.7653, 1118.4824, 1118.5532, 1118.7997, 1121.9397, 1138.5283, 1145.0836, 1154.1222, 1224.0252, 1280.9892, 1286.3355, 1293.7174, 1304.3843, 1304.4249, 1307.1626, 1307.7894, 1333.8135, 1352.5493, 1402.936, 1463.1459, 2886.2576, 2897.014, 2897.5548, 2898.0773, 2904.9758, 2906.6594, 3022.3193, 3022.3517, 3029.3245, 3029.3492, 3037.506, 3037.5529],
     'mulliken_charges': OneIndexedArray([-0.271682, 0.090648, 0.090012, 0.090649, -0.18851, 0.095355, 0.09536, -0.200782, 0.098551, 0.098567, -0.18851, 0.095364, 0.095351, -0.271682, 0.090649, 0.090012,  0.090649])
     }

Therefore, we are taking the last geometry and molecular properties from each file
and combining them into a ``ConformationalEnsemble``.

Now, let's extract out just the filename and energies using standard slicing syntax::

    property_names = ["filename", "energy"]
    conformational_energies = conformational_ensemble[:,property_names]

We can then determine the lowest energy and display the results in a `pandas` dataframe::

    df = DataFrame(conformational_energies, columns=property_names)
    df["rel_energy"] = (df.energy - df.energy.min()) * 627.509469
    print(df)

The output is::

                                     filename    energy  rel_energy
    0  test/static/pentane_conformation_1.out -0.055241    0.000000
    1  test/static/pentane_conformation_2.out -0.054881    0.226124
    2  test/static/pentane_conformation_3.out -0.054171    0.671446
    3  test/static/pentane_conformation_4.out -0.053083    1.354009

That's it!  You can find this code as a unit test (``test/test_pentane.py``).  For further
recipes and documentation, please read on!

