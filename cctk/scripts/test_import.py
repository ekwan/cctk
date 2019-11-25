import sys
import os
import numpy as np

sys.path.append(os.path.relpath('../cctk'))

from cctk import MOL2File

file = MOL2File.read_file("cctk/scripts/dodecane.mol2")

print(file.molecules.molecules[0].atomic_numbers)
