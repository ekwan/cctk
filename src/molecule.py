import sys
import re
import numpy as np

from geom_utility import distance, angle, dihedral

class Molecule():
    '''
    id = int, for use with conformational ensembles, optional
    name = string, for identification, optional
    theory = dict, containing information about theory, to keep track for basis set scans, optional
    atoms = list of atomic symbols
    coordinates  = list of 3-tuples of xyz coordinates
    bonds = Graph object containing connectivity information
    energy = holds the calculated energy of the molecule, optional
    '''
    def __init__(self, atoms, coordinates, name=None, theory=None, id=None, bonds=None):
        pass

    def formula():
        pass

    def write_gaussian_input(header, footer=None):
        pass
    
    def append_to_gaussian_input(header, footer=None):
        '''
        Appends to GaussianInputFile object; to be appended in .gjf file using Link1 command
        '''
        pass

    def write_mol2_input(header, footer=None):
        pass

    def set_distance():
        pass

    def set_angle():
        pass

    def set_dihedral():
        pass

    def calculate_mass_spectrum():
        pass
    
    def calculate_bonds_from_vdw():
        pass

    def add_atom_at_centroid(symbol, atom_numbers):
        pass
