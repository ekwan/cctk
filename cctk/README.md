## cctk (main directory)

## Contents: 

### Core Classes: 

`molecule.py`: methods for dealing with molecules, including changing structure and assigning connectivity

`helper_functions.py`: provides helper methods for calculating dihedral angles, getting atomic numbers, and so forth

`ensemble.py`: methods for dealing with groups of molecules

`input_file.py` and `output_file.py`: generic classes for input and output files

### Gaussian I/O: 

`gaussian_output_file.py`: represents a single Gaussian output file

`gaussian_data.py`: factory methods to create a Gaussian output file object

`parse_gaussian.py`: helper methods to extract key terms from Gaussian `.out` files

`gaussian_input_file.py`: represents a single Gaussian input file

`gaussian_job.py`: factory methods to create a Gaussian input file object
