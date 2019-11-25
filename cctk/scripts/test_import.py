import sys
import os
import numpy as np

sys.path.append(os.path.relpath('../cctk'))

from cctk import MOL2File

file = MOL2File.read_file("cctk/scripts/dodecane.mol2")

mol = file.molecules[0]

print(mol.bonds.edges())

mol.assign_connectivity()

print(mol.bonds.edges())
