.. _recipe_09:

======================================
Reading and Writing ORCA Files
======================================

- ``import cctk`` is assumed.
- *cctk* was originally designed with Gaussian in mind, but supports basic writing of ORCA input files and parsing of ORCA output
- The default print level in ORCA (`! NormalPrint`) is recommended for parsing with cctk
- Note that `! MiniPrint` may make Orca output files unreadable by cctk
- All of the input and output files used in these recipes can be accessed `here <./../test/static/14-butanedione.gjf>`_.
- All of the input and output files used in these recipes can be accessed `there <https://github.com/ekwan/cctk/tree/master/test/static>`_.
- All of input files used in these recipes can also be copied as text at the end of this recipe.

"""""""""""""""""""""""""""""""""""""""
Writing a simple ORCA input file
"""""""""""""""""""""""""""""""""""""""


- In this recipe, we convert an ``.xyz`` file into an ORCA ``.inp`` file.
- cctk can only write single-geometry input files for ORCA

::

    read_path = "test/static/test_peptide.xyz"
    write_path = "test/static/test_peptide.inp"

    file = cctk.XYZFile.read_file(read_path)
    header = "! aug-cc-pVTZ aug-cc-pVTZ/C DLPNO-CCSD(T) TightSCF TightPNO MiniPrint\n%pal nproc 4 end\n%maxcore 4000\n%mdci\n    density none\nend"
    cctk.OrcaFile.write_molecule_to_file(write_path, file.molecule, header)

This returns the following ``.inp`` file.

::

  ! aug-cc-pVTZ aug-cc-pVTZ/C DLPNO-CCSD(T) TightSCF TightPNO
  %maxcore 4000
  %pal
    nproc 4
  end
  %mdci
    density none
  end

  * xyz 0 1
  7         -2.59192634   -2.21198726    0.00000000
  1         -1.59236634   -2.24165320    0.00000000
  6         -3.27888036   -0.93617529    0.00000000
  1         -3.90371323   -0.85937226   -0.88982302
  6         -4.16550350   -0.79005325    1.23214102
  6         -2.28840137    0.21943574    0.00000000
  1         -3.55273032   -0.85570127    2.13119602
  8         -1.05147529   -0.01204927    0.00000000
  7         -2.77667761    1.60597324    0.00000000
  1         -3.75940657    1.78999817    0.01945589
  6         -1.83400333    2.70605183   -0.02819160
  1         -1.18843234    2.65457511    0.84855705
  6         -0.96291101    2.65067720   -1.27869356
  6         -2.55550909    4.04615879   -0.02290161
  1         -1.59452200    2.71041346   -2.16503406
  8         -3.81307602    4.08358955    0.00347381
  8         -1.80434930    5.26272535   -0.04790345
  6         -2.67382050    6.36615276    0.21924435
  1         -2.10181522    7.26956749    0.25863054
  1         -3.40495658    6.44054794   -0.55844808
  1         -3.16527271    6.21398306    1.15744436
  9         -0.27272248    1.49057508   -1.29610038
  9         -0.09860228    3.68770504   -1.27210140
  9         -4.78867817    0.40706316    1.19949698
  9         -5.08726215   -1.77632380    1.24357510
  6         -3.36436987   -3.46267939    0.00000000
  8         -4.58360243   -3.47888422   -0.13970006
  6         -2.54201651   -4.75124693    0.18689264
  1         -2.19710803   -5.09555006   -0.76567924
  1         -1.70233405   -4.54997635    0.81881291
  1         -3.15523267   -5.50351715    0.63739824
  *

""""""""""""""""""""""""""""""""""""""""""""""""
Writing an ORCA input file from a SMILES string
""""""""""""""""""""""""""""""""""""""""""""""""

- In this recipe we use the SMILES string of uridine to quickly generate an input file for geometry optimization and frequency calculation.
- Writing molecules from SMILES requires RDKIT which can be installed with ``pip install rdkit``

::

    import cctk, sys

    write_path = "test/static/orca_uridine_opt_freq.inp"

    # define a smiles string
    SMILES = "C1=CN(C(=O)NC1=O)C2C(C(C(O2)CO)O)O"

    # define a cctk molecule from SMILES string
    mol = cctk.Molecule.new_from_smiles(SMILES)

    cctk.OrcaFile.write_molecule_to_file(write_path, mol, 
	    header="! b3lyp/G 6-31g(d) D3 CPCM(water) opt freq tightscf",
	    variables={"maxcore": 1000},
	    blocks={"pal": ["nproc 4"] 
            },
        )

This returns the following ``.inp`` file.

::

    ! b3lyp/G 6-31g(d) D3 CPCM(water) opt freq tightscf
    %maxcore 1000
    %pal
      nproc 4
    end

    * xyz 0 1
    6          3.81472230    0.33750522    0.14614943
    6          2.70351934    0.51809502   -0.57053703
    7          1.44507718    0.26401237   -0.05548295
    6          1.27698672   -0.21741368    1.23580658
    8          0.18983726   -0.48697883    1.74334991
    7          2.42703629   -0.39430982    1.96055412
    6          3.69690514   -0.15269515    1.52668941
    8          4.68861151   -0.33199659    2.22232127
    6          0.27652907    0.46311659   -0.91630781
    6         -0.43325704   -0.82726794   -1.36208439
    6         -1.75885904   -0.78542173   -0.61961108
    6         -1.98108745    0.71387446   -0.44165945
    8         -0.67931598    1.29470444   -0.23261185
    6         -2.88835597    1.07023466    0.73521751
    8         -2.50963593    0.39274159    1.93051863
    8         -2.78649902   -1.41911137   -1.38349736
    8         -0.67447364   -0.80085635   -2.78132772
    1          4.80783558    0.53296161   -0.23537906
    1          2.76381254    0.87833512   -1.59442019
    1          2.32244110   -0.73365402    2.90389299
    1          0.58795124    1.03504050   -1.79944336
    1          0.13832243   -1.73754263   -1.15744293
    1         -1.65814829   -1.31425822    0.33219174
    1         -2.39890432    1.15967953   -1.35317397
    1         -3.92426777    0.79417652    0.51661742
    1         -2.84833598    2.14471841    0.94008678
    1         -1.52917957    0.41925484    1.99037051
    1         -3.54789329   -1.53307223   -0.78323817
    1         -1.52137470   -1.28387237   -2.89754915
    *

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Reading Simple ORCA Output Files
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Reading an ORCA output file is as simple as::

  path = 'test/static/orca_uridine_opt_freq.out'
  file = cctk.OrcaFile.read_file(path)

From the resulting ORCA file we can acccess several properties that apply to the whole job::
  
  file.job_types

returns a list of job types::

  [<OrcaJobType.OPT: 'opt'>, <OrcaJobType.FREQ: 'freq'>, <OrcaJobType.SP: 'sp'>]

indicating the job contains an opt and freq calculations. A single point calculation is included in all jobs::

  file.successful_terminations
  # returns '3' confirming successful termination of all three job types

  file.num_imaginaries()
  # returns '0', confirming that there are no imaginary frequencies in the optimzied structure

We can also access the full ensemble of geometries from the geometry optimization and their corresponding properties::

  path = 'test/static/orca_uridine_opt_freq.out'
  file = cctk.OrcaFile.read_file(path)

  file.ensemble.properties_list()
  # returns a list of dictionaries
  # each dictionary in the list corresponds to a geometry from the optimization
  # each dictionary contains property keys mapped to property values for the specified geometry.

  file.ensemble.get_properties_list()[0]

Returns the properties dictionary of the first geometry in the ensemble::

  # {'energy': -911.078699594417,
  # 'filename': './test/static/orca_uridine_opt_freq.out',
  # 'iteration': 0,
  # 'scf_iterations': 13.0,
  # 'rms_gradient': 0.004391566,
  # 'max_gradient': 0.0237232031,
  # 'rms_step': 0.0198286671,
  # 'max_step': 0.067347534}

To access the frequencies from the last geometry in the ensemble::

  freqs = file.ensemble[-1,'frequencies']
  # assigns a list of frequencies to the variable freqs

To access a given property for each member of the ensemble::

  geom_iters = file.ensemble[:,'iteration']
  energy = file.ensemble[:, 'energy']
  rms_grad = file.ensemble[:, 'rms_gradient']


We can then easily plot the property as a function of optimization step:: 

  import matplotlib.pyplot as plt

  f1 = plt.figure(figsize=(8,6))
  plt.scatter(geom_iters, energy)
  plt.ylabel(f"energy (hartree)")
  plt.xlabel(f"geometry step")
  plt.close()

  f2 = plt.figure(figsize=(8,6))
  plt.scatter(geom_iters, rms_gradient)
  plt.ylabel(f"rms_gradient")
  plt.xlabel(f"geometry step")
  plt.close()

Calling ``f1`` returns:

.. image:: ./img/r09_step_vs_energy.png
    :width: 450
    :align: center

Calling ``f2`` returns:

.. image:: ./img/r09_step_vs_rms_grad.png
    :width: 450
    :align: center


""""""""""""""""""""""""""""""""
Sample input files
""""""""""""""""""""""""""""""""

test_peptide.xyz::

  31
  test_peptide.xyz
  N         -2.59192634   -2.21198726    0.00000000
  H         -1.59236634   -2.24165320    0.00000000
  C         -3.27888036   -0.93617529    0.00000000
  H         -3.90371323   -0.85937226   -0.88982302
  C         -4.16550350   -0.79005325    1.23214102
  C         -2.28840137    0.21943574    0.00000000
  H         -3.55273032   -0.85570127    2.13119602
  O         -1.05147529   -0.01204927    0.00000000
  N         -2.77667761    1.60597324    0.00000000
  H         -3.75940657    1.78999817    0.01945589
  C         -1.83400333    2.70605183   -0.02819160
  H         -1.18843234    2.65457511    0.84855705
  C         -0.96291101    2.65067720   -1.27869356
  C         -2.55550909    4.04615879   -0.02290161
  H         -1.59452200    2.71041346   -2.16503406
  O         -3.81307602    4.08358955    0.00347381
  O         -1.80434930    5.26272535   -0.04790345
  C         -2.67382050    6.36615276    0.21924435
  H         -2.10181522    7.26956749    0.25863054
  H         -3.40495658    6.44054794   -0.55844808
  H         -3.16527271    6.21398306    1.15744436
  F         -0.27272248    1.49057508   -1.29610038
  F         -0.09860228    3.68770504   -1.27210140
  F         -4.78867817    0.40706316    1.19949698
  F         -5.08726215   -1.77632380    1.24357510
  C         -3.36436987   -3.46267939    0.00000000
  O         -4.58360243   -3.47888422   -0.13970006
  C         -2.54201651   -4.75124693    0.18689264
  H         -2.19710803   -5.09555006   -0.76567924
  H         -1.70233405   -4.54997635    0.81881291
  H         -3.15523267   -5.50351715    0.63739824

orca_uridine_opt_freq.inp::

  ! b3lyp/G 6-31g(d) D3 CPCM(water) opt freq tightscf
  %maxcore 1000
  %pal
    nproc 4
  end

  * xyz 0 1
  6          3.81472230    0.33750522    0.14614943
  6          2.70351934    0.51809502   -0.57053703
  7          1.44507718    0.26401237   -0.05548295
  6          1.27698672   -0.21741368    1.23580658
  8          0.18983726   -0.48697883    1.74334991
  7          2.42703629   -0.39430982    1.96055412
  6          3.69690514   -0.15269515    1.52668941
  8          4.68861151   -0.33199659    2.22232127
  6          0.27652907    0.46311659   -0.91630781
  6         -0.43325704   -0.82726794   -1.36208439
  6         -1.75885904   -0.78542173   -0.61961108
  6         -1.98108745    0.71387446   -0.44165945
  8         -0.67931598    1.29470444   -0.23261185
  6         -2.88835597    1.07023466    0.73521751
  8         -2.50963593    0.39274159    1.93051863
  8         -2.78649902   -1.41911137   -1.38349736
  8         -0.67447364   -0.80085635   -2.78132772
  1          4.80783558    0.53296161   -0.23537906
  1          2.76381254    0.87833512   -1.59442019
  1          2.32244110   -0.73365402    2.90389299
  1          0.58795124    1.03504050   -1.79944336
  1          0.13832243   -1.73754263   -1.15744293
  1         -1.65814829   -1.31425822    0.33219174
  1         -2.39890432    1.15967953   -1.35317397
  1         -3.92426777    0.79417652    0.51661742
  1         -2.84833598    2.14471841    0.94008678
  1         -1.52917957    0.41925484    1.99037051
  1         -3.54789329   -1.53307223   -0.78323817
  1         -1.52137470   -1.28387237   -2.89754915
  *


