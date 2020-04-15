.. _recipe_01:

==================================
Reading and Writing Gaussian Files
==================================

Here is how to read and write Gaussian files.  ``import cctk`` is assumed.
Statements like ``file.title = "title"`` indicate what you would see if you
printed the fields; they should not be taken literally.

"""""""""""""""""""""""""""""
Reading a Gaussian Input File
"""""""""""""""""""""""""""""

::

    # read the input file
    path = "test/static/gaussian_file.gjf"
    file = cctk.GaussianFile.read_file(path)

    # what's in the file object
    file.route_card = "#p opt freq=noraman m062x/6-31g(d) scrf=(smd,solvent=diethylether)"
    file.job_types = [cctk.JobType.OPT, cctk.JobType.FREQ, cctk.JobType.SP]
    file.link0 = {"mem": "1GB", "chk": "test.chk"}
    file.title = "title"
    file.footer = None

    # get the input geometry
    molecule = file.get_molecule()   # returns the last (and only) molecule

    # what's in the molecule object
    assert isinstance(molecule, cctk.molecule.Molecule)
    assert molecule.num_atoms() == 31
    mol.charge = 0
    mol.multiplicity = 1

::

""""""""""""""""""""""""""""""
Reading a Gaussian Output File
""""""""""""""""""""""""""""""

::

    # read the output file
    path = "test/static/gaussian_file.out"
    file = cctk.GaussianFile.read_file(path)

    # what's in the file object
    file.route_card = "#p opt freq=noraman m062x/6-31g(d) scrf=(smd,solvent=diethylether)"
    file.link0 = {"mem": "32GB",  "nprocshared": "16"}
    file.job_types = [cctk.JobType.OPT, cctk.JobType.FREQ, cctk.JobType.SP]
    file.title = "title"
    file.footer = None
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



"""""""""""""""""""""""
Multiple Link1 Sections
"""""""""""""""""""""""
    with open(old_path) as old:
        with open(new_path) as new:
            self.assertListEqual(
                list(new),
                list(old)
            )

    os.remove(new_path)

    def test_link1_out_file(self):
        path = "test/static/ethane.out"
        f, lines = cctk.GaussianFile.read_file(path, return_lines=True)

        self.assertEqual(len(lines), 3)
        self.assertEqual(len(f), 3)
        self.assertTrue(all(isinstance(file, cctk.GaussianFile) for file in f))

        self.assertListEqual(f[0].job_types, [cctk.JobType.OPT, cctk.JobType.FREQ, cctk.JobType.SP])
        self.assertListEqual(f[1].job_types, [cctk.JobType.NMR, cctk.JobType.SP])
        self.assertListEqual(f[2].job_types, [cctk.JobType.NMR, cctk.JobType.SP])


    def test_write_ensemble(self):
        path = "test/static/gaussian_file.out"
        file = cctk.GaussianFile.read_file(path)
        ense = file.ensemble

        old_path = "test/static/ensemble.gjf"
        new_path = "test/static/new_ensemble.gjf"
        cctk.GaussianFile.write_ensemble_to_file(new_path, ense, "#p opt freq=noraman b3lyp/6-31g(d)", print_symbol=True)

        with open(old_path) as old:
            with open(new_path) as new:
                self.assertListEqual(
                    list(new),
                    list(old)
                )

        os.remove(new_path)

    def test_force_extraction(self):
        path = "test/static/dcm_force.out"
        file = cctk.GaussianFile.read_file(path)
        ense = file.ensemble

        self.assertListEqual(list(ense[0, "forces"][1]), [2.672010074,2.672010074,0.0])

    def test_charges(self):
        path = "test/static/dcm_force.out"
        file = cctk.GaussianFile.read_file(path)
        ense = file.ensemble
        self.assertEqual(ense[-1, "mulliken_charges"][1], -0.051271)

        path = "test/static/h2o.out"
        file = cctk.GaussianFile.read_file(path)
        ense = file.ensemble
        self.assertEqual(ense[-1, "hirshfeld_charges"][1], -0.312885)
        self.assertEqual(ense[-1, "hirshfeld_spins"][1], 0)

    def test_dipole(self):
        path = "test/static/dcm_force.out"
        file = cctk.GaussianFile.read_file(path)
        ense = file.ensemble
        self.assertEqual(ense[-1, "dipole_moment"], 0.3316)

