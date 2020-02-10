.. _ensemble:

============
Ensembles
============

An ensemble in *cctk* is simply a collection of molecules. If desired, an ensemble can have energies associated with it. 

Molecules can be added to an ensemble through ``add_molecule()``, which uses Python's ``copy.deepcopy()`` method to create a copy of the molecule being added (to prevent unexpected side-effects).

Molecules in an ensemble can be accessed like a list::

    ensemble = cctk.Ensemble()

    ensemble.add_molecule(mol1)
    ensemble.add_molecule(mol2)

    ensemble[0]     #### returns mol1
    
    ensemble[0] = mol2 
    
    ensemble[0]     #### returns mol2

Conformational Ensembles
=========================

Conformational ensemble objects inherit from ensembles, so they can be manipulated in the same way. 
Unlike ensembles, however, conformational ensembles impose the restriction that every molecule must have the same atomic numbers in the same order. 
Therefore, a conformational ensemble can be used to describe:

- Points along the reaction coordinate
- Different conformations
- Different spatial arrangements of several molecules

Conformational ensembles can be "aligned" to orient certain atoms in the same direction using the ``align()`` method, and redundant conformations can be removed with ``eliminated_redundant()`` 
(which simply compares the RMS deviation between each geometry and deletes those below a threshold). 

To see the lowest energy structures, the ``get_lowest_energy()`` or ``get_within_cutoff()`` methods can be used. 
