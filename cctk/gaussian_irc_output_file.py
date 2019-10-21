import sys
import re
import numpy as np

class GaussianIRCOutputFile(GaussianOutputFile):
    '''
    energies = list of energies for each cycle
    scf_iterations = number of iterations per cycle
    geometries = list of geometrries for each cycle
    done = Boolean
    '''
    def __init__(self):
        pass    

    def read_file(filename,text):
        super.read_file()

    def create_molecule_from_geometry(number):
        pass 

    def print_energies():
        pass

    def print_geometric_parameters():
        for geometry in self.geometries():
            pass
