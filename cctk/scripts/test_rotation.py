import sys
import os

sys.path.append(os.path.relpath('../cctk'))

from cctk import Molecule
from cctk.helper_functions import to_degrees, to_radians

#H3 = Molecule(atoms=[20,20,20,1], geometry=[[1,0,0],[0,0,0],[0,1,0],[4,4,4]])
H3 = Molecule(atoms=[1,1,1,1], geometry=[[-1,-1,1],[0,0,1],[0,1,1],[0,0,0]])
#H3 = Molecule(atoms=[20,20,20], geometry=[[1,0,0],[0,0,0],[0,1,0]])
print(H3.bonds.edges)
H3.set_angle(1,2,3,20)
print(H3.geometry)
