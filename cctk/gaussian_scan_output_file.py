import sys
import re
import numpy as np

from cctk import GaussianOutputFile

class GaussianScanOutputFile(GaussianOutputFile):
    """
    Class for Gaussian output files originating from a geometry optimization. 
    
    Attributes:
        atoms (list): list of atoms 
        energies (list): list of energies for each cycle
        scf_iterations (list): number of iterations per cycle
        max_displacements (list): list of max displacement values for each cycle 
        rms_displacements (list): list of rms displacement values for each cycle 
        max_forces (list): list of max force values for each cycle 
        rms_forces (list): list of rms force values for each cycle 
        gradient (list): list of gradient values for each cycle 
        geometries (list): list of geometries for each cycle (so really a list of lists)
    """
    def __init__(self):
        pass    

    def read_file(filename,text):
        super.read_file()

    def get_lowest_gradient_geometry():
        pass
        
    def create_molecule_at_lowest_gradient():
        mol_geom = self.get_lowest_gradient_geometry()
        pass

    def create_molecule_from_geometry(number):
        pass 

    def print_energies():
        pass

    def print_geometric_parameters():
        for geometry in self.geometries():
            pass
