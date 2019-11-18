from collections import OrderedDict
import numpy as np
from mol2 import read_mol2
from mae import read_mae

# read a mae (macromodel) file with conformers in it
geometries, symbols, bonds, property_names, property_values, contains_conformers = read_mae("dodecane_csearch-out.mae")
geometry_index = -1
for symbol,position in zip(symbols[geometry_index],geometries[geometry_index]):
    x,y,z = position
    print(f"{symbol:2s} {x:12.6f} {y:12.6f} {z:12.6f}")

property_dic_zero = OrderedDict(zip(property_names[geometry_index], property_values[geometry_index]))
for property_name,property_value in property_dic_zero.items():
    print(f"{property_name:50s} : {property_value:50s}")

print(bonds[geometry_index].edges)
#print(contains_conformers)

exit()
input("press return to continue")

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


