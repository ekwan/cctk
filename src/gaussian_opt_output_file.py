import sys
import re
import numpy as np

class GaussianOptOutputFile(GaussianOutputFile):
    '''
    energies = list of energies for each cycle
    scf_iterations = number of iterations per cycle
    rms_forces
    max_forces
    rms_displacements
    max_displacements
    geometries = list of geometrries for each cycle
    done = Boolean
    '''
    def __init__(self):
        pass    

    def read_file(filename,text):
        super.read_file()

    def get_lowest_energy_geometry():
        pass
        
    def create_molecule()
        mol_geom = self.get_lowest_energy_geometry()
        pass

    def create_molecule_from_geometry(number):
        pass 

    def print_energies():
        pass

    def print_geometric_parameters():
        for geometry in self.geometries():
            pass
