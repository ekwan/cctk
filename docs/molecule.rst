.. _molecule:

=========
Molecules
=========

The ``Molecule`` object is central to *cctk*, since every input and output file contains one. 

Molecules in *cctk* have two principal attributes: ``atomic_numbers`` (a list of atomic numbers) and ``geometry`` (a list of Cartesian coordinates).
Both attributes are stored internally as ``numpy`` arrays, but input lists will automatically be cast to ``numpy``::

    atom_nums = [1, 1]
    geometry = [[0, 0, 0], [1, 0, 0]]

    dihydrogen = cctk.Molecule(atom_nums, geometry)

Molecules can also be given charge (``charge``) and spin multiplicity (``multiplicity``); these attributes are not directly manipulated, but may be required to generate certain filetypes. 

Indexing
========

The indexing of atoms is problematic since Python inherently 0-indexes its data structures, but most computational chemistry software uses 1-indexing. 
To prevent confusion when switching between various programs, *cctk* has adopted the convention of 1-indexing atom numbers for all outward-facing methods. 
Internal methods (prefixed with an underscore) are sometimes 0-indexed and sometimes 1-indexed; if you suspect something is amiss, check the relevant documentation. 

Wherever possible, 1-indexed accessor methods have been added to preclude the need to directly interface with 0-indexed Python code::

    molecule.geometry[4]        #### intrinsic Python indexing; will return the values for atom 5
    molecule.atomic_numbers[4]

    molecule.get_vector(4)      #### accessor methods; will return the values for atom 4
    molecule.get_atomic_number(4)

Bonding and Connectivity
========================

*cctk* stores bond information as a ``networkx`` graph, with bond order stored as the weight of the graph edge. 

If the bond connectivity is known, then a simple list of bonds can be read in to ``Molecule.__init__()``::

    atom_nums = [54, 9, 9]
    geometry = [[0, 0, 0], [0, 0, 1], [0, 0, -1]]
    bonds = [[1, 2], [1, 3]]

    xef2 = cctk.Molecule(atom_nums, geometry, bonds)

If, however, the bond connectivity is not known (e.g. in ``.xyz`` files), then *cctk* can try to predict the connectivity using ``assign_connectivity()``.
This simply compares the distance between every set of atoms to a table of covalent radii, so it's not guaranteed to work, but for "normal" organic molecules it has a high success rate. 

Once created, the bonding graph can be queried with ``get_bond_order(atom1, atom2)``. 
To get a list of an atom's neighbors, call ``get_adjacent_atoms(atom)``; to see if two atoms are in the same molecule, try ``are_connected(atom1, atom2)``::

    xef2.get_adjacent_atoms(1)  #### returns [2, 3]
    xef2.get_adjacent_atoms(2)  #### returns [1]

    xef2.are_connected(2, 3)    #### returns True 

If you want to explicitly change the bonding, you can use the ``add_bond()`` or ``remove_bond()`` methods. 

Working with the Molecule
=========================

Individual atoms can be added or removed with ``add_atom()`` or ``remove_atom()`` methods, and atoms can be added more precisely using ``add_atom_at_centroid()``.

To access precise geometric properties, one can use ``get_distance()``, ``get_angle()``, and ``get_dihedral()`` (the corresponding ``set_`` methods all modify these properties). 
Specific dihedral angles can be "optimized" (by minimizing the pairwise atomic RMSD) using the relatively slow method ``optimize_dihedral()``.

To manipulate the entire molecule, one can use ``rotate_molecule()`` and ``translate_molecule()``.

Finally, every atom in the molecule can be "wiggled" through addition of random noise using the ``perturb()`` method: this can be useful for avoiding spurious minima or collinearity. 
