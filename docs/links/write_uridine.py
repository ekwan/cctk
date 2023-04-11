import cctk, sys

write_path = "../../test/static/orca_uridine_opt_freq.inp"

# define a smiles string
SMILES = "C1=CN(C(=O)NC1=O)C2C(C(C(O2)CO)O)O"

# define a cctk molecule from SMILES string
mol = cctk.Molecule.new_from_smiles(SMILES)

cctk.OrcaFile.write_molecule_to_file(write_path, mol, 
    header="! b3lyp/G 6-31g(d) D3 CPCM(water) opt freq tightscf")