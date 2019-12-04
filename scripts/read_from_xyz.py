import sys
from cctk import XYZFile, GaussianFile

#### This file reads an XYZ file and outputs a Gaussian input file. 

#### Usage: ``python read_from_xyz.py my_molecule`` (to read in ``my_molecule.xyz``)

filename = sys.argv[1]

file = XYZFile.read_file(f"{filename}.xyz")
file.molecule.assign_connectivity()

GaussianFile.write_molecule_to_file(
    f"{filename}.gjf", 
    file.molecule, 
    "#p opt=(ts,calcfc,noeigentest) freq=noraman bp86/def2tzvp empiricaldispersion=gd3bj",
    None,
)
