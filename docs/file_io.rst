.. _file_io:

=============================
File I/O
=============================

All files inherit from the *cctk* ``File`` method. 

To read in a file, use the ``read_file()`` method; to write to a file, use the ``write_file()`` method. 


Gaussian Files
==============

Reading in a Gaussian ``.out`` file using ``GaussianFile.read_file()`` creates a ``GaussianFile`` object.
The core of this object is a ``ConformationalEnsemble`` object (``GaussianFile.molecules``) which describes the output structure or structures. 

*cctk* also automatically extracts the header and stores it as ``GaussianFile.header()``.
Footer autodetection is challenging due to the wide variety of Gaussian footers, but ``opt=modredundant`` footers are autodetected and stored in ``GaussianFile.footer``.

Each output file is assigned one or more ``JobType`` instances, which are automatically detected from the header. 
Some methods are restricted to certain job types -- for instance, calling ``file.num_imaginaries()`` on a file that doesn't have type ``FREQ`` will result in an error. 

*cctk* extracts the following additional parameters from Gaussian ``.out`` files:

- The title of the file (``GaussianFile.title``)
- The number of successful terminations (``GaussianFile.success``)
- The energies at each cycle (``GaussianFile.energies``)
- The number of SCF iterations required for each cycle (``GaussianFile.scf_iterations``)
- The RMS displacement after each cycle (``GaussianFile.rms_displacements``)
- The RMS force after each cycle (``GaussianFile.rms_forces``)
- The list of frequencies, for ``FREQ`` jobs (``GaussianFile.frequencies``)
- The enthalpy, for ``FREQ`` jobs (``GaussianFile.enthalpy``)
- The Gibbs free energy, for ``FREQ`` jobs (``GaussianFile.gibbs_free_energy``)

Calling ``read_file()`` on a ``.gjf`` file will result in the creation of a Gaussian file object with only minimal information (header, footer, title, and the molecule). 

To extract a specific ``Molecule`` object from a Gaussian file, use the ``get_molecule()`` method: this returns the last molecule by default, but can be passed a specific number. 

Calling ``write_file()`` on an existing Gaussian file object will result in the creation of a ``.gjf`` file with the same header, footer, etc. 
By default the last molecule from ``GaussianFile.molecules()`` will be written, but this can be overridden.
