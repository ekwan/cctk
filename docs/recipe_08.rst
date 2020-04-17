.. _recipe_08:

==========================
Creating and Adding Groups
==========================

- ``import cctk`` is assumed.

"""""""""""""""""""""""
Using a Preloaded Group
"""""""""""""""""""""""

- *cctk* comes with several groups defined by default (listed below).
- To add a group, use the ``add_group_to_molecule`` function to generate a new molecule.
- The atom the group is added to must only be making one bond; usually a hydrogen is a good choice.

**Preloaded Groups**

=================================   ===========================================
Name                                Synonyms 
=================================   ===========================================
``methyl``                          ``Me`` ``CH3``
``ethyl``                           ``Et``, ``C2H5``
``isopropyl``                       ``iPr``, ``iC3H7``
``tert-butyl``                      ``tBu``, ``tC4H9``
``hydroxy``                         ``OH``
``methoxy``                         ``MeO``, ``OMe``, ``CH3O``
``acetamido``                       ``NHAc``
``amino``                           ``NH2``
``dimethylamino``                   ``Me2N``, ``NMe2``
``trifluoromethyl``                 ``CF3``
``cyano``                           ``CN``
``nitro``                           ``NO2``
``carboxylmethyl``                  ``MeO2C``, ``CO2Me``
``fluoro``                          ``F``
``chloro``                          ``Cl``
``bromo``                           ``Br``
``iodo``                            ``I``
``pentafluorosulfanyl``             ``SF5``
``sulfonyl``                        ``SO3H``
``acetyl``                          ``Ac``, ``COMe``
``formyl``                          ``CHO``
=================================   ===========================================

::

    from cctk.load_groups import load_group

    # load our starting molecule
    formic_acid = cctk.GaussianFile.read_file("formic_acid.out").get_molecule()

    # to load a predefined group, use the group's name
    trifluoromethyl = load_group("trifluoromethyl")
    assert isinstance(trifluoromethyl, cctk.Group)

    # now we can create a new molecule by adding our group to an existing molecule 
    # here we specify to add the group to atom 4, the formyl proton, to give us TFA
    trifluoroacetic_acid = cctk.Group.add_group_to_molecule(formic_acid, trifluoromethyl, 4)
    assert isinstance(trifluoroacetic_acid, cctk.Molecule)


""""""""""""""""""""""
Aligning Two Molecules
""""""""""""""""""""""

- We can create a new ``Group`` object from any molecule by specifying which atom to make the attachment point from. 
- As before, the ``attach_to`` atom must only be making one bond; usually a hydrogen is a good choice.

::

    # define a new group, where the attachment point is atom 6 (a methyl-group hydrogen)
    file = cctk.GaussianFile.read_file("acetaldehyde.out")
    acetaldehyde = file.get_molecule()
    group = cctk.Group.new_from_molecule(attach_to=6, molecule=acetaldehyde)

    # now create a new molecule by adding the new group to atom 5 (another hydrogen on acetaldehyde)
    new_mol = cctk.Group.add_group_to_molecule(acetaldehyde, group, 5)
    file.write_file("1,4-butanedione.gjf", molecule=new_mol)


