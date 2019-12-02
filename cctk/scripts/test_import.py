import sys
import os
import numpy as np

sys.path.append(os.path.relpath('../cctk'))

from cctk import MOL2File, MAEFile

file = MOL2File.read_file("cctk/scripts/dodecane_csearch.mol2")
#file, names, vals = MAEFile.read_file("cctk/scripts/dodecane_csearch.mae")

#energies = [x[-7].astype(float) for x in vals]
dihedrals = [None] * 8
for i in range(0, 8):
    dihedrals[i] =  file.molecules.get_geometric_parameters("dihedral", i+1, i+2, i+3,i+4)

print("Num:          Energies:     Dihedrals:")
for i, mol in enumerate(file.molecules.molecules):
    current_dihedrals = [x[i] for x in dihedrals]
    dihedral_string = "   ".join(f"{x:6.2f}" for x in current_dihedrals)
#    print(f"{i+1:6}        {energies[i]:10.2f}    {dihedral_string}")
    print(f"{i+1:6}        {dihedral_string}")

#### these energies seem bizarre...

