.. _recipe_06:

================
NMR Spectroscopy
================

"""""""""""""""
Chemical Shifts
"""""""""""""""

- ``import cctk`` is assumed.
- Isotropic shieldings are read automatically.
- The ``scale_nmr_shifts`` function will apply scaling factors to the shieldings
  to produce chemical shift predictions.
- The default scaling factors are very rough:

::

    # dict: atomic symbol --> (slope, intercept)
    # defines the slope to be positive
    DEFAULT_NMR_SCALING_FACTORS = {
            "H" : (1.0716,  31.6660),
            "C" : (1.0300, 180.4300),
            "N" : (0.9776, 244.5626)
    }

- You may provide your own scaling factor dictionary.
- Only elements for which scalings are provided will be considered.
- The ``symmetrical_atom_numbers`` parameter tells *cctk* which nuclei are
  equivalent (e.g. methyl group protons).
- Note that if symmetrical atom numbers are provided for some atoms, but
  no scaling factors are given, an error will result.

::

    # read file
    gaussian_file = cctk.GaussianFile.read_file("test/static/LSD_custom.out")
    ensemble = gaussian_file.ensemble

    # the raw shieldings
    # note: 1-indexed
    shieldings = ensemble[:,"isotropic_shielding"]

    # can automatically detect symmetric atoms from methyl, isopropyl, or tert-butyl groups
    molecule = ensemble.molecules[-1].assign_connectivity()
    symmetric = molecule.get_symmetric_atoms() # [[32, 33, 34], [37, 38, 39]]

    # linearly scale to get shifts
    scaled_shifts, shift_labels = cctk.helper_functions.scale_nmr_shifts(
        ensemble,
        symmetrical_atom_numbers=symmetric, 
        scaling_factors="default"
    )

    expected_shifts = [6.52352,6.6285,6.51045,6.53005,6.22303,2.11021,2.7025,2.73022,2.38541,2.35172,3.1467,5.82979,
                       3.29202,1.92326,114.924,98.3836,107.641,94.3333,104.421,109.795,95.1041,112.168,121.346,
                       45.4898,14.1014,26.7028,36.3779,29.4323,104.708,155.804,38.0661,109.579,22.7099]
    expected_shifts = np.asarray(expected_shifts)
    assert np.abs(scaled_shifts[0] - expected_shifts) <= 0.001).all()

    # shift labels are provided to facilitate downstream analysis
    # order parallels that of scaled_shifts
    # note: list of lists!
	shift_labels == [['H20' 'H21' 'H22' 'H23' 'H24' 'H25' 'H26' 'H27' 'H28' 'H29' 'H30' 'H31'
	                  'H37/38/39' 'H32/33/34' 'C1' 'C2' 'C3' 'C4' 'C6' 'C7' 'C8' 'C9' 'C10'
	                  'C11' 'C12' 'C14' 'C15' 'C16' 'C17' 'C18' 'C36' 'N5' 'N13']]

""""""""""""""""""
Coupling Constants
""""""""""""""""""

- The final J couplings from an ``nmr=spinspin`` or ``nmr=mixed`` calculation will be automatically parsed.
- Data are stored in a ``j_couplings`` in ``properties_dict``.
- This is a symmetric 2D ``np.array`` where the two axes represent the number of atoms.  Diagonal elements
  are 0.  Each value is given in Hz.
- As a result, the couplings array is **zero-indexed in both dimensions**.
- Note that ``nmr=mixed`` occurs in two internal job steps in the same ``Link1``.
- In all cases of coupling constant calculations, the ``j_couplings`` property can be found in the *last*
  ``properties_dict``.

::

    # this is a single point nmr=mixed calculation
    # as a result, there is one ``Link`` section, but there are two (identical) geometries
    gaussian_file = cctk.GaussianFile.read_file("test/static/acetone-couplings1.out")
    ensemble = gaussian_file.ensemble

    # couplings and shieldings will be found in the second (and last) dictionary
    shieldings = ensemble[-1,"isotropic_shielding"]
    expected_shieldings = [165.8515, 30.794, 30.93, 30.9302, -21.4514, -375.1462,
                           159.4249, 30.6991, 30.6993, 30.8559]
    self.assertTrue((np.abs(shieldings - expected_shieldings) <= 0.0001).all())

    # couplings are provided as a 2D np.array
    expected_couplings = np.array(\
    [[  0. ,124.1,134.7,134.7, 34.8, -0.8, 15.5, -0.4, -0.4,  4.9],
     [124.1,  0. ,-14.4,-14.4, -3.7, -2.1,  0.6,  0.4,  0.4,  0.7],
     [134.7,-14.4,  0. ,-20.4, -6.9, -1.3,  1.1, -0.6, -1.3, -0.1],
     [134.7,-14.4,-20.4,  0. , -6.9, -1.3,  1.1, -1.2, -0.6, -0.1],
     [ 34.8, -3.7, -6.9, -6.9,  0. , 43.7, 35. , -5.7, -5.6, -6.4],
     [ -0.8, -2.1, -1.3, -1.3, 43.7,  0. , -1.1, -1.9, -1.9, -0.8],
     [ 15.5,  0.6,  1.1,  1.1, 35. , -1.1,  0. ,127.2,127.2,137.5],
     [ -0.4,  0.4, -0.6, -1.2, -5.7, -1.9,127.2,  0. ,-19.1,-14.6],
     [ -0.4,  0.4, -1.3, -0.6, -5.6, -1.9,127.2,-19.1,  0. ,-14.6],
     [  4.9,  0.7, -0.1, -0.1, -6.4, -0.8,137.5,-14.6,-14.6,  0. ]])
    couplings = ensemble[-1,"j_couplings"]
    self.assertTrue(np.any(expected_couplings-couplings < 0.1))



