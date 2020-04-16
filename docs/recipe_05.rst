.. _recipe_05:

==========
NMR Shifts
==========

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

::

    # read file
    gaussian_file = cctk.GaussianFile.read_file("test/static/LSD_custom.out")
    ensemble = gaussian_file.ensemble

    # the raw shieldings
    # note: 1-indexed
    shieldings = ensemble[:,"isotropic_shielding"]

    # linearly scale to get shifts
    scaled_shifts, shift_labels = cctk.helper_functions.scale_nmr_shifts(ensemble,
                                  symmetrical_atom_numbers=[[37,38,39],[32,33,34]], scaling_factors="default")
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
