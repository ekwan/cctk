.. _recipe_01:

==================================
Reading and Writing Gaussian Files
==================================

- ``import cctk`` is assumed.
- Statements like ``file.title == "title"`` or ``assert molecule.num_atoms() == 31``
  indicate what you would see if you printed the fields.

"""""""""""""""""""""""""""""
Reading a Gaussian Input File
"""""""""""""""""""""""""""""

- ``file.get_molecule()`` is equivalent to ``file.ensemble.molecules[-1]``.

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

- **Important: only files specifying verbose output with** ``#p`` **in the route card
  will be parsed correctly.**

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

"""""""""""""""""""""""""""""""""""""""""""""
Writing One Molecule to a Gaussian Input File
"""""""""""""""""""""""""""""""""""""""""""""

- Only route cards specifying verbose output (``#p``) are allowed to
  ensure compatibility with *cctk*.

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

- Each link will be parsed into a separate `GaussianFile`.
- The `properties_dict` key `link1_idx1` identifies which ``Link1`` the geometry came from.

::

    # read a file with multiple link1 directives
    path = "test/static/ethane.out"
    files = cctk.GaussianFile.read_file(path)

    # get back a list of file objects
    len(files) == 3
    for file in files:
        assert isinstance(file, cctk.GaussianFile)

    # different links can correspond to different types of jobs
    files[0].job_types == [cctk.JobType.OPT, cctk.JobType.FREQ, cctk.JobType.SP]
    files[1].job_types == [cctk.JobType.NMR, cctk.JobType.SP]
    files[2].job_types == [cctk.JobType.NMR, cctk.JobType.SP]

"""""""""""""""""""""""""""""""""""""""""""""""""""""
Writing Multiple Molecules to One Gaussian Input File
"""""""""""""""""""""""""""""""""""""""""""""""""""""

- The geometries will be combined into a single job using ``Link1``.
- Here, a single route card is specified.  However, list-like inputs can be specified
  for more complex jobs.  See the API documentation for details.
- The ``footer`` specifies what goes after each geometry.  This can be useful for specifying
  special basis sets.
- The ``print_symbol`` flag specifies whether elements should be specified by atomic
  number or symbol.

::

    assert isinstance(ensemble, cctk.Ensemble)
    cctk.GaussianFile.write_ensemble_to_file(filename, ensemble, route_card = "#p opt freq=noraman b3lyp/6-31g(d)",
                                             link0={"mem": "32GB", "nprocshared": 16}, footer=None,
                                             title="title", print_symbol=False)

""""""""""""""""""""""""""""""""""""""""""""""""""""
Using Custom Basis Sets from the Basis Set Exchange
""""""""""""""""""""""""""""""""""""""""""""""""""""

- Bespoke basis sets can be downloaded automatically from the `Basis Set Exchange <https://www.basissetexchange.org/>`.
- By default, the ``add_custom_basis_set`` method appends the basis set to the footer. However, 
  passing the ``return_string`` option allows for increased control over formatting (e.g. for combination with ``opt=modredundant``).
- The ``gen`` keyword should be used in combination with these basis sets.

::

    assert isinstance(file, cctk.GaussianFile)
    file.route_card = "#p opt wB97XD/gen"
    file.add_custom_basis_set("pcseg-2")

    assert isinstance(file2, cctk.GaussianFile)
    file2.route_card = "#p opt=modredundant wB97XD/gen"
    basis = file2.add_custom_basis_set("pcseg-2", return_string=True)
    file2.footer = f"B 1 10 F\n\n{basis}"

