import sys
import re
import numpy as np

from geom_utility import distance, angle, dihedral, align_molecules

class Ensemble():
    '''
    name = string, for identification
    molecules = list of molecules
    '''

    def __init__(self, name=None): 
        pass    

    def eliminate_redundant():
        pass    

    def lowest_energy(self, n=1):
        pass    
    
    def write_gaussian_input():
        for molecule in molecules:
            molecule.write_gjf_input()    

    def write_mol2_input():
        for molecule in molecules:
            molecule.write_mol2_input()    

    def set_distance():
        for molecule in molecules:
            molecule.set_distance()    

    def set_angle():
        for molecule in molecules:
            molecule.set_angle()    

    def set_dihedral():
        for molecule in molecules:
            molecule.set_dihedral()    
