import cctk, sys

# this is a script for generating gjf files from molecule names
# usage: gen_start_orca.py mol_name file_name.inp
# example: gen_start_orca.py FK506 fk506.inp

mol = cctk.Molecule.new_from_name(sys.argv[1])
cctk.OrcaFile.write_molecule_to_file(sys.argv[2], mol, 
	header="! opt b3lyp/G 6-31g(d) D3 Normalprint Printbasis PrintMOs #CPCM",
	variables={"maxcore": 1000},
	blocks={"pal": ["nproc 4"] 
            # , "mdci": ["density none"]
            # , "cpcm": ["smd true", "SMDsolvent \"dichloromethane\""]
            # , "scf": ["Print[P_SCFMemInfo] 1"]
            },
        )
print("done")
