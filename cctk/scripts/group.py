import sys
import os
import numpy as np
import copy

sys.path.append(os.path.relpath('../cctk'))

from cctk import GaussianData, Molecule, Group

output_file = GaussianData.read_opt('cctk/scripts/acetaldehyde.out')

molecule = Molecule(output_file.atoms, output_file.get_final_geometry())
molecule.assign_connectivity()

print(molecule.geometry)

group = Group.new_from_molecule(attach_to=6, molecule=molecule)
print(group.geometry)

