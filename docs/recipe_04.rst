.. _recipe_04:

=============================
Bonds and Molecular Fragments
=============================

- ``import cctk`` is assumed.

"""""""""""""""""
Bond Connectivity
"""""""""""""""""

- `cctk` keeps track of bonded atoms in a `networkx <https://https://networkx.github.io/>`_ graph.
- Bond orders are not tracked (yet).
- If connectivity information is not available, it can be automatically assigned based on atomic radii.

::

    # read a Molecule
    path="test/static/test_peptide.xyz"
    molecule = cctk.XYZFile.read_file(path).molecule
    
    # the atomic numbers
    molecule.atomic_numbers == [7, 1, 6, 1, 6, 6, 1, 8, 7, 1, 6, 1, 6, 6, 1, 8, 8, 6, 1, 1, 1, 9, 9, 9, 9, 6, 8, 6, 1, 1, 1]

    # xyz files have no connectivity, so automatically assign
    molecule.assign_connectivity()
    molecule.bonds.edges() = [(1, 2), (1, 3), (1, 26), (3, 4), (3, 5), (3, 6), (5, 7), (5, 24), (5, 25), (6, 8), (6, 9), (9, 10), (9, 11), (11, 12), (11, 13), (11, 14), (13, 15), (13, 22), (13, 23), (14, 16), (14, 17), (17, 18), (18, 19), (18, 20), (18, 21), (26, 27), (26, 28), (28, 29), (28, 30), (28, 31)]


"""""""""""""""""""
Molecular Fragments
"""""""""""""""""""

::

        (frag1, frag2) = mol._get_bond_fragments(3, 5)
        self.assertEqual(len(frag1), 27)
        self.assertEqual(len(frag2), 4)

        self.assertEqual(len(mol._get_fragment_containing(5)), 31)
        mol.remove_bond(3,5)
        self.assertEqual(len(mol._get_fragment_containing(5)), 4)
        self.assertFalse(mol.are_connected(3,5))
        mol.add_bond(3,5)
        self.assertEqual(len(mol._get_fragment_containing(5)), 31)
        self.assertTrue(mol.are_connected(3,5))

""""""""""""
Adding Atoms
""""""""""""

::

    def test_add_atoms(self):
        mol = cctk.Molecule(np.array([2], dtype=np.int8), [[0, 0, 0]])
        self.assertEqual(mol.num_atoms(), 1)

        mol.add_atom("He", [1, 0, 0])
        self.assertListEqual(mol.atomic_numbers.tolist(), [2, 2])
        self.assertEqual(mol.num_atoms(), 2)

        mol.add_atom("Ar", [3, 0, 0])
        self.assertEqual(mol.num_atoms(), 3)

        mol.add_atom_at_centroid("He", [2, 3])
        self.assertEqual(mol.num_atoms(), 4)
        self.assertListEqual(list(mol.get_vector(4)), [2, 0, 0])

""""""""""""""
Atom Selection
""""""""""""""

::

    def test_selection(self):
        mol = self.load_molecule()
        self.assertListEqual(mol.get_heavy_atoms(), [1, 3, 5, 6, 8, 9, 11, 13, 14, 16, 17, 18, 22, 23, 24, 25, 26, 27, 28])
        self.assertListEqual(mol.get_atoms_by_symbol("F"), [22, 23, 24, 25])

"""""""""""""""""""""
Translating Molecules
""""""""""""""""""""""

::

        mol = cctk.Molecule(np.array([12], dtype=np.int8), [[0, 0, 0]])

        v = np.array([1.5234,1.231234,-1.77777])
        mol = mol.translate_molecule(v)

        self.assertListEqual(mol.geometry.tolist()[0], list(v))
        self.assertTrue(isinstance(mol.geometry, cctk.OneIndexedArray))

        mol2 = self.load_molecule()
        v2 = np.zeros(shape=3)

        mol2_shift = mol2.translate_molecule(v2)
        self.assertListEqual(mol2.geometry.tolist()[0], mol2_shift.geometry.tolist()[0])

"""""""""""""""""""
Combining Molecules
"""""""""""""""""""

::

        m1 = cctk.Molecule(np.array([12], dtype=np.int8), [[0, 0, 0]], charge=-1, multiplicity=1)
        m2 = cctk.Molecule(np.array([12], dtype=np.int8), [[2, 0, 0]], charge=2, multiplicity=1)

        m3 = cctk.Molecule.combine_molecules(m1, m2)
        self.assertTrue(isinstance(m3, cctk.Molecule))
        self.assertEqual(m3.num_atoms(), 2)
        self.assertEqual(m3.charge, 1)
        self.assertEqual(m3.multiplicity, 1)


""""""""""""""""
Molecular Volume
""""""""""""""""

::

        mol = self.load_molecule()
        self.assertEqual(mol.volume(), 80.42662712363737)
    
""""""""""""
Mass Spectra
""""""""""""

::

        mol = cctk.Molecule(np.array([12], dtype=np.int8), [[0, 0, 0]])
        masses, weights = mol.calculate_mass_spectrum()
        self.assertListEqual(list(masses), [23.])
