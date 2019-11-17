import numpy as np
from mol2 import read_mol2
#import mae

# read a single mol2
geometry, symbols, bonds, contains_conformers = read_mol2("dodecane.mol2")
print(np.shape(geometry))
print(geometry)
print(np.shape(symbols))
print(symbols)
print(np.shape(bonds))
print(bonds.edges)
print(contains_conformers)
input("press return to continue")

# read a mol2 with conformers in it
geometries, symbols, bonds, contains_conformers = read_mol2("dodecane-csearch.mol2", contains_conformers=True)
print(np.shape(geometries))
print(geometries[0])
print(geometries[-1])
print(np.shape(symbols))
print(symbols)
print(np.shape(bonds))
print(bonds.edges)
print(contains_conformers)
input("press return to continue")


# read a mae (macromodel) file with conformers in it
