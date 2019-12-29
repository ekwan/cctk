import sys
import os
import numpy as np
import copy

sys.path.append(os.path.relpath('../cctk'))

from cctk import GaussianFile, Molecule, Group, XYZFile
from cctk.group_substitution import add_group_to_molecule

output_file = GaussianFile.read_file('cctk/scripts/acetaldehyde.out')

#molecule.assign_connectivity()

print(output_file.get_molecule().geometry)

group = Group.new_from_molecule(attach_to=6, molecule=output_file.get_molecule())
print(group.geometry)

new_mol = add_group_to_molecule(output_file.get_molecule(), group, 5)
print(new_mol.geometry)

output_file.write_file("cctk/scripts/1,4-butanedione.gjf", molecule=new_mol)
