import sys
import os
import numpy as np

sys.path.append(os.path.relpath('../cctk'))

from cctk import XYZData, Molecule, Ensemble, GaussianJob

# write a bunch of output files
output_file = XYZData.read_xyz('cctk/scripts/CpG.xyz')

molecule = Molecule(output_file.atoms, output_file.geometry)
molecule.assign_connectivity()
ensemble = Ensemble(name='cpg conformers')

angles = [0, 60, 120, 180, 240, 240.5]
for x in angles:
    for y in angles:
        molecule.set_dihedral(1, 7, 6, 8, x)
        molecule.set_dihedral(23, 24, 25, 1, y)

        ensemble.add_molecule(molecule)

#ensemble.align(atoms=[1, 53, 25])
ensemble.eliminate_redundant()

for molecule in ensemble.molecules:
    molecule.check_for_conflicts() 
    input_file = GaussianJob.create_opt(molecule.atoms, molecule.geometry)
    x = molecule.get_dihedral(1, 7, 6, 8)
    y = molecule.get_dihedral(23, 24, 25, 1)
    input_file.write_file(filename=f"cctk/scripts/CpG_conformers/CpG_{int(round(x))}_{int(round(y))}.gjf")

#output_file.write_file('cctk/scripts/CpG_copy.xyz')
