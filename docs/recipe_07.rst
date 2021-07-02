.. _recipe_07:

======================================
Reading and Writing Other File Formats
======================================

- ``import cctk`` is assumed.
- *cctk* was originally designed with Gaussian in mind, but has some limited
  support for other file formats.  Please consider `OpenBabel <http://openbabel.org/wiki/Main_Page>`_
  for more extensive file conversion needs.

"""""""""
XYZ Files
"""""""""

- *cctk* assumes the first line is the number of atoms, the second line is the title,
  and the third and subsequent lines are geometry specifications of the form
  ``atom_symbol x_position y_position z_position``.
- For xyz files containing multiple structures, concatenate blocks of the above form
  without any separating blank lines.

::

    # read xyz
    path = "test/static/test_peptide.xyz"
    file = cctk.XYZFile.read_file(path)
    
    # file contents
    file.titles[0] == "peptide example"
    ensemble = file.ensemble

    # convenient accessor
    molecule = file.get_molecule()

    # write xyz
    assert isinstance(molecule2, Molecule)
    XYZFile.write_molecule_to_file(filename, molecule2, title="title"):

""""""""""
MOL2 Files
""""""""""

- Both single and multiple structures per ``.mol2`` are supported.

::

    # read a single structure
    path = "test/static/dodecane.mol2"
    file = cctk.MOL2File.read_file(path)
    assert isinstance(file, cctk.MOL2File)
    ensemble = file.ensemble
    assert isinstance(ensemble, cctk.ConformationalEnsemble)
    len(ensemble) == 1
    molecule = ensemble.molecules[0]

    # read multiple structures
    path = "test/static/dodecane-csearch.mol2"
    file = cctk.MOL2File.read_file(path)
    assert isinstance(file, cctk.MOL2File)
    ensemble = file.ensemble
    assert isinstance(ensemble, cctk.ConformationalEnsemble)
    len(ensemble) == 597

    # write one structure
    filename2 = "my_filename.mol2"
    cctk.MOL2File.write_molecule_to_file(filename2, molecule, title)

    # write multiple structures
    filename3 = "your_filename.mol2"
    cctk.MOL2File.write_ensemble_to_file(filename3, ensemble)

- Setting ``print_status_messages`` to ``True`` will print progress to stdout.
- For large files that contain only conformers, set ``contains_conformers`` to ``True``.
- This will skip many consistency checks and avoid the redundant copying of many
  data structures.  This speeds up performance by over 10x.

::

    # faster reading of MOL2 files with conformers
    file = cctk.MOL2File.read_file(path, print_status_messages=False, contains_conformers=True)

"""""""""""""
Maestro Files
"""""""""""""

- Only reading is supported.
- Both single and multiple structures are supported.

::

    # read .mae file
    path = "test/static/dodecane_csearch-out.mae"
    file, pnames, pvals = cctk.MAEFile.read_file(path)
    assert isinstance(file, cctk.MAEFile)
   
    # file contents
    ensemble = file.ensemble
    assert isinstance(ensemble, cctk.ConformationalEnsemble)
    len(ensemble) == 597
    len(pnames) == 597
    len(pvals) == 597

""""""""""""""""
Orca Input Files
""""""""""""""""

- In this recipe, we convert an ``.xyz`` file into an ``.inp`` file for Orca.
- Only single geometries are supported.

::

    read_path = "test/static/test_peptide.xyz"
    write_path = "test/static/test_peptide.inp"

    file = cctk.XYZFile.read_file(read_path)
    header = "! aug-cc-pVTZ aug-cc-pVTZ/C DLPNO-CCSD(T) TightSCF TightPNO MiniPrint\n%pal nproc 4 end\n%maxcore 4000\n%mdci\n    density none\nend"
    cctk.OrcaFile.write_molecule_to_file(write_path, file.molecule, header)



