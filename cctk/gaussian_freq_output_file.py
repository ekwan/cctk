import sys
import re
import numpy as np

from cctk import GaussianOutputFile

class GaussianFreqOutputFile(GaussianOutputFile):
    '''
    Class for Gaussian output files originating from a geometry optimization. 
    
    Attributes:
        atoms (list): list of atoms 
        energies (list): list of energies for each cycle
        scf_iterations (list): number of iterations per cycle
        max_displacements (list): list of max displacement values for each cycle 
        rms_displacements (list): list of rms displacement values for each cycle 
        max_forces (list): list of max force values for each cycle 
        rms_forces (list): list of rms force values for each cycle 
        geometries (list): list of geometries for each cycle (so really a list of lists)
        frequencies (list): list of frequencies
        gibbs_free_energy (float): gibbs free energy, from vibrational correction
        enthalpy (float): enthalpy, from vibrational correction
    '''
    def __init__(self):
        pass    

    def read_file(filename,text):
        super.read_file()

    def create_molecule():
        mol_geom = self.geometry()
        pass

    def num_imaginary():
        pass
