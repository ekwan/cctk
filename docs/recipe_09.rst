.. _recipe_09:

======================================
Reading and Writing ORCA Files
======================================

- ``import cctk`` is assumed.
- *cctk* was originally designed with Gaussian in mind, but has some limited
  support for other file formats.  Please consider `OpenBabel <http://openbabel.org/wiki/Main_Page>`_
  for more extensive file conversion needs.


""""""""""""""""
Writing a simple ORCA input file
""""""""""""""""

- In this recipe, we convert an ``.xyz`` file into an ``.inp`` file for Orca.
- Only single geometries are supported.
- Note that the MiniPrint option may make Orca output files unreadable by cctk

::

    read_path = "test/static/test_peptide.xyz"
    write_path = "test/static/test_peptide.inp"

    file = cctk.XYZFile.read_file(read_path)
    header = "! aug-cc-pVTZ aug-cc-pVTZ/C DLPNO-CCSD(T) TightSCF TightPNO MiniPrint\n%pal nproc 4 end\n%maxcore 4000\n%mdci\n    density none\nend"
    cctk.OrcaFile.write_molecule_to_file(write_path, file.molecule, header)

""""""""""""""""
Writing an ORCA input file with multiple jobs from one molecule
""""""""""""""""


""""""""""""""""
Writing a compound job to find a stationary point
""""""""""""""""


""""""""""""""""
Writing a compound job to find a transition state with one imaginary frequency
""""""""""""""""


""""""""""""""""
Writing Multiple Molecules to One Gaussian Input File
""""""""""""""""


""""""""""""""""
Creating Input Files by Molecule Name or SMILES
""""""""""""""""
