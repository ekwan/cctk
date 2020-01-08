import sys
from cctk import XYZFile, GaussianFile

filename = sys.argv[1]

file = XYZFile.read_file(f"{filename}.xyz")

file.molecule.assign_connectivity()

GaussianFile.write_molecule_to_file(
    f"{filename}.gjf",
    file.molecule,
    "#p opt=(ts,calcfc,noeigentest) freq=noraman bp86/def2tzvp empiricaldispersion=gd3bj",
    None,
)
