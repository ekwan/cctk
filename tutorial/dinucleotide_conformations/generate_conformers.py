from cctk import XYZFile, ConformationalEnsemble, GaussianFile

#### This file generates a bunch of different phosphate-based rotamers for the given CpG dinucleotide.
#### Usage: ``python generate_conformers.py``

#### The input file is in the ".xyz" format, which means we have to infer the bonds ourselves (unlike other file formats).
output_file = XYZFile.read_file('CpG.xyz')
output_file.molecule.assign_connectivity()

#### Here we create a blank ``ConformationalEnsemble`` object to hold the new conformers.
ensemble = ConformationalEnsemble(name='cpg conformers')

#### For each of the defined rotatable bonds, we're going to set them to every value in ``angles``. 
angles = [0, 60, 120, 180, 240, 240.4, 360]
for x in angles:
    for y in angles:
        output_file.molecule.set_dihedral(1, 7, 6, 8, x)
        output_file.molecule.set_dihedral(23, 24, 25, 1, y)
        ensemble.add_molecule(output_file.molecule)

#### ``ensemble.eliminate_redundant()`` will not only catch the fact that 0 and 360 are the same angle, but also notice that 240 and 240.4 aren't very different. 
#### The cutoff can be adjusted higher or lower depending on the situation. 

old_num = len(ensemble.molecules)
ensemble.eliminate_redundant()
new_num = len(ensemble.molecules)
print(f"originally {old_num} conformers, but after eliminating redundant there are {new_num}!")

#### Now all that's left to do is to take the resultant molecules and write them to an output file, with the appropriate naming scheme. 
for molecule in ensemble.molecules:
    molecule.check_for_conflicts()
    x = molecule.get_dihedral(1, 7, 6, 8)
    y = molecule.get_dihedral(23, 24, 25, 1)
    GaussianFile.write_molecule_to_file(f"conformers/CpG_{int(round(x))}_{int(round(y))}.gjf", molecule, "#p opt b3lyp/6-31g(d)", None)
