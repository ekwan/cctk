from cctk import XYZFile, ConformationalEnsemble, GaussianFile

#### This file generates a bunch of different phosphate-based rotamers for the given CpG dinucleotide.
#### Usage: ``python generate_conformers.py``

output_file = XYZFile.read_file("CpG.xyz")
mol = output_file.molecule.assign_connectivity()

ensemble = ConformationalEnsemble()

angles = [0, 60, 120, 180, 240, 241, 300]
for x in angles:
    for y in angles:
        mol.set_dihedral(1, 7, 6, 8, x)
        mol.set_dihedral(23, 24, 25, 1, y)
        ensemble.add_molecule(mol, copy=True)

old_num = len(ensemble)
ensemble = ensemble.eliminate_redundant()
new_num = len(ensemble)
print(f"originally {old_num} conformers, but after eliminating redundant there are {new_num}!")

count = 0
for molecule in ensemble.molecules:
    x = int(round(molecule.get_dihedral(1, 7, 6, 8)))
    y = int(round(molecule.get_dihedral(23, 24, 25, 1)))
    if molecule.check_for_conflicts():
        GaussianFile.write_molecule_to_file(f"CpG_{x}_{y}.gjf", molecule, route_card="#p opt b3lyp/6-31g(d)")
        count += 1
    else:
        print(f"x={x} y={y} no good - error")

print(f"wrote {count} molecules to files")
