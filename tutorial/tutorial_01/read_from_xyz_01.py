import re
from cctk import XYZFile, GaussianFile, OrcaFile

filename = "./tutorial1.xyz"
file = XYZFile.read_file(filename)
newfile = filename.rsplit('/',1)[-1]
newfile = re.sub(r"xyz$", "gjf", newfile)

molecule = file.get_molecule()

GaussianFile.write_molecule_to_file(
    newfile,
    molecule,
    "#p opt freq=noraman b3lyp/6-31g(d) empiricaldispersion=gd3bj",
)

# to write and orca input simultaneously we could use the block below
newfile = re.sub(r"gjf$", "inp", newfile)
OrcaFile.write_molecule_to_file(
    newfile,
    molecule,
    "! opt freq b3lyp/6-31g(d) d3bj",
    )