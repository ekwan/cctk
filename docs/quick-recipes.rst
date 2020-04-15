.. _quick-recipes:

=============
Quick Recipes
=============

------------------------
Conformations of Pentane
------------------------

*How to determine the lowest energy structure from a set of Gaussian output files.*

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

