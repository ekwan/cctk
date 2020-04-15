.. _recipe_01:

==================================
Reading and Writing Gaussian Files
==================================

Here is how to read and write Gaussian files.  ``import cctk`` is assumed.
Statements like ``file.title == "title"`` indicate what you would see if you
printed the fields.

"""""""""""""""""""""""""""""
Reading a Gaussian Input File
"""""""""""""""""""""""""""""

::

    # read the input file
    path = "test/static/gaussian_file.gjf"
    file = cctk.GaussianFile.read_file(path)

    # what's in the file object
    file.route_card == "#p opt freq=noraman m062x/6-31g(d) scrf=(smd,solvent=diethylether)"
    file.job_types == [cctk.JobType.OPT, cctk.JobType.FREQ, cctk.JobType.SP]
    file.link0 == {"mem": "1GB", "chk": "test.chk"}
    file.title == "title"
    file.footer == None

    # get the input geometry
    molecule == file.get_molecule()   # returns the last (and only) molecule

    # what's in the molecule object
    assert isinstance(molecule, cctk.molecule.Molecule)
    assert molecule.num_atoms() == 31
    mol.charge == 0
    mol.multiplicity == 1

""""""""""""""""""""""""""""""
Reading a Gaussian Output File
""""""""""""""""""""""""""""""

::

    # read the output file
    path = "test/static/gaussian_file.out"
    file = cctk.GaussianFile.read_file(path)

    # what's in the file object
    file.route_card == "#p opt freq=noraman m062x/6-31g(d) scrf=(smd,solvent=diethylether)"
    file.link0 == {"mem": "32GB",  "nprocshared": "16"}
    file.job_types == [cctk.JobType.OPT, cctk.JobType.FREQ, cctk.JobType.SP]
    file.title == "title"
    file.footer == None
    assert isinstance(file.ensemble, cctk.ConformationalEnsemble)

    # following is equivalent to file.get_molecule()
    ensemble = file.ensemble
    molecule = file.molecules[-1]

    # get the molecular properties dictionary
    properties_dict = ensemble.get_properties_dict(molecule)
    properties_dict["filename"] == path
    properties_dict["energy"] == -1159.56782622

"""""""""""""""""""""""""""""
Writing a Gaussian Input File
"""""""""""""""""""""""""""""

::

    # need an initial molecule object
    assert isinstance(molecule, cctk.Molecule)

    # define the options for the new file
    new_path = "input_file.gjf"
    link0 = {"chk": "checkpoint.chk", "mem": "32GB", "nprocshared": "16"}
    route_card = "#p opt=(ts,calcfc,noeigentest) m062x/6-31g(d)"

    # write the file
    cctk.GaussianFile.write_molecule_to_file(new_path, molecule, route_card, link0)


"""""""""""""""""""""""
Multiple Link1 Sections
"""""""""""""""""""""""

::

    # read a file with multiple link1 directives
    path = "test/static/ethane.out"
    files = cctk.GaussianFile.read_file(path)

    # get back a list of file objects
    self.assertEqual(len(files), 3)
    self.assertTrue(all(isinstance(file, cctk.GaussianFile) for file in f))

    # different files have different types and route_cards
    self.assertListEqual(f[0].job_types, [cctk.JobType.OPT, cctk.JobType.FREQ, cctk.JobType.SP])
    self.assertListEqual(f[1].job_types, [cctk.JobType.NMR, cctk.JobType.SP])
    self.assertListEqual(f[2].job_types, [cctk.JobType.NMR, cctk.JobType.SP])

    # we can also write multiple molecules to the same input file using link1
    assert isinstance(ensemble, cctk.Ensemble)
    cctk.GaussianFile.write_ensemble_to_file(new_path, ensemble, "#p opt freq=noraman b3lyp/6-31g(d)")

