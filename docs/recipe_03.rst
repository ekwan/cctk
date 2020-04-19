.. _recipe_03:

==========================================
Measuring and Setting Internal Coordinates
==========================================

- ``import cctk`` is assumed
- Assume a ``Molecule`` has been loaded as ``molecule``.
- *cctk* assumes molecular geometries are given in Angstroms

""""""""""""""""""""""""""""""
Measuring Distances and Angles
""""""""""""""""""""""""""""""

::

    # the atomic numbers as a one-indexed array
    atomic_numbers = molecule.atomic_numbers

    # measure bond distance between atoms 1 and 2 in A
    bond_distance = molecule.get_distance(1,2)

    # measure bond angle between atoms 1, 2, and 3 in degrees
    bond_angle = molecule.get_angle(1, 2, 3)

    # measure dihedral angle between atoms 1, 2, 3, and 4 in degrees
    dihedral_angle = molecule.get_dihedral(1, 2, 3, 4)

    # alternately, all get_ or set_ methods will also take a list of atom numbers
    dihedral_angle = molecule.get_dihedral(atoms=[1, 2, 3, 4])

""""""""""""""""""""""""""""
Setting Distances and Angles
""""""""""""""""""""""""""""

- ``Molecule`` objects are modified in-place.
- When ``move`` is set to a ``group``, rings are not allowed (see below).

::

    # make a copy of the molecule if desired
    molecule2 = copy.deepcopy(molecule)

    # set bond distance between atoms 1 and 2 to 2.00 A
    # move="group": all atoms connected to atom 2 will also be moved (this is the default)
    #               if atom 2 is connected to any atom that is also connected to atom1,
    #               there is a ring and an error will be thrown
    # move="atom": only atom 2 will be moved
    molecule2.set_distance(1, 2, 2.00, move="group")
    
    # set the bond angle between atoms 1, 3, and 5 to 120 degrees
    # move="group": all atoms connected to the third atom will also be moved (default, rings not allowed)
    # move="atom": only the third atom will be moved
    molecule2.set_angle(1, 3, 5, 120)

    # set the dihedral angle between atoms 1, 7, 2, and 8 to 157 degrees
    # move="atom": only the fourth atom will be moved
    # move="group4": all atoms bonded to the fourth atom will be moved (rings not allowed)
    # move="group34": all atoms bonded to the third and fourth atoms will be moved (default, rings not allowed)
    molecule2.set_dihedral(1, 7, 2, 8, 157)

    # use unpacking to pass in the atom numbers as a tuple or list
    atom_tuple = (1, 7, 2, 8)
    molecule2.set_dihedral(*atom_tuple, 157)
    atom_list = [1, 7, 2, 8]
    molecule2.set_dihedral(*atom_list, 157)

    # or just pass the list directly!
    molecule2.set_dihedral(atoms=atom_list)

""""""""""""""""""""
Checking For Clashes
""""""""""""""""""""

- Checks for steric clashes based on covalent radii.
- If two atoms are closer than the sum of their vdw radii + ``min_buffer``, they are said to be clashing.
- Returns ``True`` if there are no conflicts, or ``False`` if there are.

::

        no_clashes_present = molecule.check_for_conflicts()
        

