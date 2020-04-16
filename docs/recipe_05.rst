.. _recipe_05:

=======================================================
Aligning Structures and Redundant Conformer Elimination
=======================================================

- ``import cctk`` is assumed.

""""""""""""""""
Calculating RMSD
""""""""""""""""

- We can calculate the root-mean-square deviation between two conformers.
- Calculating the RMSD between two non-conformers is not supported (yet).
- No alignment is performed; the two geometries are compared as-is.

::

    # load molecule
    path = "test/static/gaussian_file.out"
    gaussian_file = cctk.GaussianFile.read_file(path)
    ensemble = gaussian_file.ensemble
    m1 = ensemble.molecules[0]
    m2 = ensemble.molecules[-1]
    
    # measure RMSD in Angstroms
    RMSD = cctk.helper_functions.compute_RMSD(m1,m2)
    RMSD == 0.0006419131435567976

""""""""""""""""""""""
Aligning Two Molecules
""""""""""""""""""""""

- We can align all the conformers in an ensemble onto one of it members.
- All conformers are translated to the origin.  Then, a different rotation is applied to each
  conformer until the RMSD between the conformer and the template conformer is minimized.
- A new ``ConformationalEnsemble`` is created and the original ``Molecule`` objects are not modified.
- The `Kabsch algorithm <https://en.wikipedia.org/wiki/Kabsch_algorithm>`_ is used.
- Alignment of non-conformers is not supported (yet).

::

    # create a ConformationalEnsemble
    path = "test/static/phenylpropane*.out"
    conformational_ensemble = cctk.ConformationalEnsemble()
    for filename in sorted(glob.glob(path)):
        gaussian_file = cctk.GaussianFile.read_file(filename)
        ensemble = gaussian_file.ensemble
        molecule = ensemble.molecules[-1]
        properties_dict = ensemble.get_properties_dict(molecule)
        conformational_ensemble.add_molecule(molecule,properties_dict)

    # define which atoms will be used for alignment
    # alternatively, set comparison_atoms to "heavy" or "all"
    comparison_atoms = [1,2,3,4,5,6]

    # perform the alignment
    # to_geometry is the template geometry (0-indexed)
    # set compute_RMSD to True to calculate the RMSD before and after the aligment
    aligned_ensemble, before_RMSD, after_RMSD = conformational_ensemble.align(to_geometry=np.int64(0), comparison_atoms=comparison_atoms, compute_RMSD=True)

    # write out the result to a Gaussian input file
    # open in GaussView and select "Single new molecule group for all files" in the "Open Files" dialog box under "Target"
    # select geometries of interest of "select all" to compare the aligned molecules
    cctk.GaussianFile.write_ensemble_to_file("test/static/phenylpropane_aligned.gjf", aligned_ensemble, "#p")

    # alternatively, write all structures to a .mol2 file and open in Maestro or a similar program
    # it be necessary to run molecule.assign_connectivity() on each member of the Ensemble first
    # if the bond connectivity information is not available
    cctk.MOL2File.write_ensemble_to_file("test/static/phenylpropane_aligned.mol2", aligned_ensemble)


"""""""""""""""""""""""""""""""
Redundant Conformer Elimination
"""""""""""""""""""""""""""""""

- Useful for removing redundant conformers after conformational searches.
- Every geometry is aligned and then a new ``ConformationalEnsemble`` is created that only contains
  unique conformations.
- Whenever redundancies are found, the lowest energy conformer is kept.
- The current ``ConformationalEnsemble`` is not modified.

::

    # need a ConformationalEnsemble
    assert isinstance(conformational_ensemble, ConformationalEnsemble)

    # eliminate redundant conformers
    # RMSD_cutoff is in Angstroms
    # comparison_atoms is a list of atom numbers, "heavy" (default), or "all"
    unique_conformers = conformational_ensemble.eliminated_redundant(RMSD_cutoff=0.5, comparison_atoms="heavy")

    # how many conformers are left
    len(unique_conformers)

    # write out the entire ensemble
    # alternatively, index the ensemble to retrieve a subset of the structures first
    cctk.GaussianFile.write_ensemble_to_file("test/static/phenylpropane_aligned2.gjf", ensemble2, "#p")
