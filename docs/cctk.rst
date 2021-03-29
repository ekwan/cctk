.. _cctk:

.. |br| raw:: html

  <br/>

.. currentmodule:: cctk

cctk package
============

This section describes the *cctk* API.  Please program responsibly.

"""""""""""
*cctk* Core
"""""""""""

*These are the core classes that represent molecules and collections of molecules in cctk*

.. autosummary::
    :nosignatures:
    :toctree: _autosummary

    Molecule <molecule.Molecule>
    Ensemble <ensemble.Ensemble>
    ConformationalEnsemble <ensemble.ConformationalEnsemble>
 
"""""
Files
"""""

*These classes enable input/output functions.*

.. autosummary::
    :nosignatures:
    :toctree: _autosummary 

    File <file.File>
    GaussianFile <gaussian_file.GaussianFile>
    MAEFile <mae_file.MAEFile>
    MOL2File <mol2_file.MOL2File>
    OrcaFile <orca_file.OrcaFile>
    PDBFile <pdb_file.PDBFile>
    SIFile <si_file.SIFile>
    XYZFile <xyz_file.XYZFile>

"""""""""""""
Miscellaneous
"""""""""""""

.. autosummary::
    :nosignatures:
    :toctree: _autosummary
    
    Group <group.Group>
    VibrationalMode <vibrational_mode.VibrationalMode>
    helper_functions
    OneIndexedArray <array.OneIndexedArray>
    LazyLineObject <lines.LazyLineObject>
    
|br|

