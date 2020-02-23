from cctk import XYZFile, ConformationalEnsemble, GaussianFile

#### This file generates a bunch of different phosphate-based rotamers for the given CpG dinucleotide.
#### Usage: ``python generate_conformers.py``

output_file = XYZFile.read_file("CpG.xyz")
output_file.molecule.assign_connectivity()

ensemble = ConformationalEnsemble(name='cpg conformers')

angles = [0, 60, 120, 180, 240, 241, 300]
for x in angles:
    for y in angles:
        output_file.molecule.set_dihedral(1, 7, 6, 8, x)
        output_file.molecule.set_dihedral(23, 24, 25, 1, y)
        ensemble.add_molecule(output_file.molecule)

old_num = len(ensemble.molecules)
ensemble = ensemble.eliminate_redundant()
new_num = len(ensemble.molecules)
print(f"originally {old_num} conformers, but after eliminating redundant there are {new_num}!")

count = 0
for molecule in ensemble.molecules:
    x = int(round(molecule.get_dihedral(1, 7, 6, 8)))
    y = int(round(molecule.get_dihedral(23, 24, 25, 1)))
    try:
        molecule.check_for_conflicts()
        GaussianFile.write_molecule_to_file(f"conformers/CpG_{x}_{y}.gjf", molecule, "#p opt b3lyp/6-31g(d)", None)
        count += 1
    except ValueError as e:
        print(f"x={x} y={y} no good - {e}")

print(f"wrote {count} molecules to files")
