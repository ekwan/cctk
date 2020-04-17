.. _recipe_04:

====================================
Bond Connectivity and Atom Selection
====================================

- ``import cctk`` is assumed.

"""""""""""""""""
Bond Connectivity
"""""""""""""""""

- `cctk` keeps track of bonded atoms in a `networkx <https://https://networkx.github.io/>`_ graph.
- This graph is stored as ``molecule.bonds``.
- In principle, bond orders are supported by calling ``bonds[atom1][atom2]["weight"]``.
  However, the existing parsers assume a bond order of 1, regardless of what is actually in the file.
- If connectivity information is not available, it can be automatically assigned based on atomic radii.

::

    # read a Molecule
    path="test/static/test_peptide.xyz"
    molecule = cctk.XYZFile.read_file(path).molecule
    
    # the atomic numbers
    molecule.atomic_numbers == [7, 1, 6, 1, 6, 6, 1, 8, 7, 1, 6, 1, 6, 6, 1, 8, 8, 6, 1, 1, 1, 9, 9, 9, 9, 6, 8, 6, 1, 1, 1]

    # xyz files have no connectivity, so automatically assign
    molecule.assign_connectivity()
    molecule.bonds.edges() == [(1, 2), (1, 3), (1, 26), (3, 4), (3, 5), (3, 6), (5, 7), (5, 24), (5, 25), (6, 8), (6, 9), (9, 10), (9, 11), (11, 12), (11, 13), (11, 14), (13, 15), (13, 22), (13, 23), (14, 16), (14, 17), (17, 18), (18, 19), (18, 20), (18, 21), (26, 27), (26, 28), (28, 29), (28, 30), (28, 31)]



""""""""""""""
Atom Selection
""""""""""""""

::

    # get atom numbers of all heavy atoms
    heavy_atom_numbers = molecule.get_heavy_atoms()
    heavy_atom_numbers == [1, 3, 5, 6, 8, 9, 11, 13, 14, 16, 17, 18, 22, 23, 24, 25, 26, 27, 28]

    # get all atom numbers for a given element
    fluorine_atom_numbers = molecule.get_atoms_by_symbol("F")
    fluorine_atom_numbers == [22, 23, 24, 25]


