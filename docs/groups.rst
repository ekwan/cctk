.. _groups:

============
Groups
============

Groups are a convenient way of combining two separate molecules into the same molecule (ordinarily a challenging task even in a 3D editor). 
To combine two molecules, one must first be converted into a ``Group``, and then ``cctk.Group.add_group_to_molecule()`` can be used to stich the two together. 


Defining Groups
===============

``Group`` objects inherit from ``Molecule`` but possess one additional property: ``attach_to``, the atom in the group which a new fragment will replace.
This must be defined upon group creation, and tells *cctk* how to combine this group with other molecules.
Internally, *cctk* also tracks the position of the atom adjacent to ``attach_to`` as ``Group.adjacent``:: 

    c2h2_atom = [6, 6, 1, 1]
    c2h2_geom = [[0, 0, 0], [0, 0, 1], [0, 0, 2], [0, 0, -1]]

    c2h2 = cctk.Molecule(c2h2_atom, c2h2_geom)
    c2h2.assign_connectivity()

    ethynyl = cctk.Group.new_from_molecule(c2h2, 3)

    ethynyl.attach_to   #### returns 3
    ethynyl.adjacent    #### returns 2 (the atom adjacent to 3)

In the above example, ``C2`` is the atom which will end up bonded to the other molecule (which will replace ``H3``).

Attempting to define an atom bonded to multiple atoms as ``attach_to`` will result in an error::

    ethynyl = cctk.Group.new_from_molecule(c2h2, 2) #### will not succeed!

Adding Groups to Molecules
==========================

Groups

Predefined Groups
=================
