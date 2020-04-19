.. _recipe_02:

======================================================================
Extracting Molecular Properties (Energies, Frequencies, Charges, etc.)
======================================================================

"""""""""""""""""""""
Example: ``opt freq``
"""""""""""""""""""""

- ``import cctk`` is assumed.
- This is a routine optimization using ``opt freq``.
- The route card job type are both parsed.

::

    # read the file
    filename = "test/static/gaussian_file.out"
    gaussian_file = cctk.GaussianFile.read_file(filename)
    
    # gaussian_file properties
    gaussian_file.route_card == "#p opt freq=noraman m062x/6-31g(d) scrf=(smd,solvent=diethylether)"
    gaussian_file.job_types == [cctk.JobType.OPT, cctk.JobType.FREQ, cctk.JobType.SP]


- An ``Ensemble`` is a collection of Molecules.
- Here, the ``Ensemble`` contains the initial, intermediate, and final
  geometries for the optimization process.

::

    # ensemble properties
    ensemble = gaussian_file.ensemble
    assert isinstance(ensemble, cctk.ConformationalEnsemble)
    len(ensemble) == 3

    # getting a molecule
    assert isinstance(ensemble.molecules[0], Molecule)

    # slicing an Ensemble returns an Ensemble
    assert isinstance(ensemble[0], cctk.ConformationalEnsemble)

"""""""""""""""""""
``properties_dict``
"""""""""""""""""""

- Each ``Molecule`` is mapped to a ``dict`` called ``properties_dict``,
  where energies and other data are stored.
- Aside: this was a simple ``opt freq`` job that contained only one ``Link1``.
  So ``read_file(filename)`` only produecs one ``GaussianFile`` and thus
  ``link1_idx`` is 0 for all geometries.

::

    # contents of the ensemble
    for molecule,properties_dict in ensemble.items():
        for property_name,value in properties_dict.items():
            print(property_name, ":", value)
        print()


- The first geometry is the input geometry

::

    '''
    energy : -1159.56782625
    scf_iterations : 13
    link1_idx : 0
    filename : test/static/gaussian_file.out
    rms_force : 9e-06
    rms_displacement : 0.000547
    '''

- The second geometry is an intermediate geometry during the optimization.
- This optimization took only one intermediate step.

::

    '''
    energy : -1159.56782622
    scf_iterations : 7
    link1_idx : 0
    filename : test/static/gaussian_file.out
    rms_force : 8e-06
    rms_displacement : 0.000489
    '''

- The third geometry is the final geometry.
- It has many more properties

::
    
    '''
    energy : -1159.56782622
    scf_iterations : 1
    link1_idx : 0
    filename : test/static/gaussian_file.out
    rms_force : 8e-06
    rms_displacement : 0.000388
    enthalpy : -1159.314817
    gibbs_free_energy : -1159.388328
    frequencies : [19.7009, 35.267, 44.7198, 49.9491, 57.0893, 64.8469, 73.0098, 81.9867, 89.8636, 98.1969, 118.548, 152.3357, 159.3717, 169.6716, 191.1281, 226.4406, 253.8972, 280.2593, 302.7976, 325.3145, 350.63, 366.305, 386.4397, 427.6294, 501.8031, 509.3823, 536.3219, 548.0141, 569.8335, 580.3367, 626.0123, 637.3627, 670.326, 687.8399, 715.4247, 812.8056, 859.8734, 916.9873, 980.1514, 1014.1109, 1023.7809, 1048.9533, 1067.319, 1076.0426, 1113.4787, 1144.4635, 1169.8286, 1181.8906, 1189.4969, 1192.3828, 1194.51, 1226.7551, 1228.0879, 1252.4422, 1274.3266, 1290.5456, 1319.7655, 1348.7007, 1353.473, 1375.5382, 1419.4037, 1428.4896, 1438.5693, 1451.7563, 1455.4477, 1486.3505, 1505.4247, 1511.1827, 1515.5803, 1518.8313, 1548.1324, 1588.3775, 1807.1221, 1834.9622, 1880.3968, 3101.1372, 3109.6277, 3116.2425, 3125.0222, 3155.476, 3177.7175, 3180.833, 3201.4795, 3209.3614, 3234.7301, 3611.0649, 3614.2603]
    mulliken_charges : [-0.655761  0.390379 -0.168647  0.247314  0.490238  0.633198  0.182115
     -0.553943 -0.670923  0.406032 -0.146986  0.254725  0.473692  0.652675
      0.202    -0.494397 -0.474991 -0.266798  0.20917   0.203205  0.200146
     -0.296932 -0.306803 -0.307007 -0.292294  0.597754 -0.534396 -0.585561
      0.214057  0.190642  0.208099]
    dipole_moment : 4.8038
    '''

- Ensemble properties can be indexed and sliced as if they were 2D arrays.
- ``None`` is used for missing data.

::

    # extracting one property value	
    ensemble[0,"energy"] == -1159.56782625

    # extracting all energies
    ensemble[:,"energy"] == [-1159.56782625, -1159.56782622, -1159.56782622]

    # enthalpy is only calculated for the final geometry of the optimization
    ensemble[:,"enthalpy"] = [None, None, -1159.314817]

"""""""""""""""""
Sorting Ensembles
"""""""""""""""""

- Ensembles can be sorted by property values (e.g., energy).
- The ordering in ``ensemble.molecules`` will be updated to reflect the new order.
- Missing entries are not allowed.
- Sorting incomparable types will result in an error.
- The result is a new ``Ensemble``.  The underlying objects are not cloned.

::

    sorted_ensemble = ensemble.sort_by("energy", ascending=False)

