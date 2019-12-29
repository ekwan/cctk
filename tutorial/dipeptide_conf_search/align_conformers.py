from cctk import XYZFile, ConformationalEnsemble, GaussianFile

output_file = XYZFile.read_file('cctk/scripts/CpG.xyz')
output_file.molecule.assign_connectivity()

ensemble = ConformationalEnsemble(name='cpg conformers')

angles = [0, 60, 120, 180, 240, 240.5]
for x in angles:
    for y in angles:
        output_file.molecule.set_dihedral(1, 7, 6, 8, x)
        output_file.molecule.set_dihedral(23, 24, 25, 1, y)
        ensemble.add_molecule(output_file.molecule)

ensemble.eliminate_redundant()

for molecule in ensemble.molecules:
    molecule.check_for_conflicts()
    x = molecule.get_dihedral(1, 7, 6, 8)
    y = molecule.get_dihedral(23, 24, 25, 1)
    GaussianFile.write_molecule_to_file(f"cctk/scripts/CpG_conformers/CpG_{int(round(x))}_{int(round(y))}.gjf", molecule, "#p opt b3lyp/6-31g(d)", None)
