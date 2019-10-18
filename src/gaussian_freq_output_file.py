import sys
import re
import numpy as np

class GaussianFreqOutputFile(GaussianOutputFile):
    '''
    scf_iterations = number of iterations
    geometry = list of geometrries for each cycle
    done = Boolean
    frequencies = list of frequencies
    gibbs_free_energy = gibbs free energy, from vibrational correction
    enthalpy = enthalpy, from vibrational correction
    '''
    def __init__(self):
        pass    

    def read_file(filename,text):
        super.read_file()

    def create_molecule()
        mol_geom = self.geometry()
        pass

    def num_imaginary():
        pass
