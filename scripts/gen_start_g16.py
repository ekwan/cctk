import cctk, sys

# this is a script for generating gjf files from molecule names
# usage: gen_start.py mol_name file_name.gjf
# example: gen_start_g16.py FK506 fk506.gjf

mol = cctk.Molecule.new_from_name(sys.argv[1])
cctk.GaussianFile.write_molecule_to_file(sys.argv[2], mol, route_card="#p opt b3lyp/6-31g(d) empiricaldispersion=gd3bj")
print("done")
