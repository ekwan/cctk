import cctk


read_path = "../../test/static/test_peptide.xyz"
path = "../../test/static/test_peptide.inp"
new_path = "test/static/test_peptide_copy.inp"
file = cctk.XYZFile.read_file(read_path)

header = "! aug-cc-pVTZ aug-cc-pVTZ/C DLPNO-CCSD(T) TightSCF TightPNO"
variables = {"maxcore": 4000}
blocks = {"pal": ["nproc 4"], "mdci": ["density none"]}

cctk.OrcaFile.write_molecule_to_file(new_path, file.get_molecule(), header, variables, blocks)
