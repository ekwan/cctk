import sys
import os
import numpy as np
import copy

sys.path.append(os.path.relpath('../cctk'))

from cctk import GaussianFile, Molecule, Group

output_file = GaussianFile.read_file('cctk/scripts/acetaldehyde.out')

#molecule.assign_connectivity()

print(output_file.get_molecule().geometry)

group = Group.new_from_molecule(attach_to=6, molecule=output_file.get_molecule())
print(group.geometry)

