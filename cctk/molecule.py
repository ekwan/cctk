import sys
import re
import numpy as np
import NetworkX as nx

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
        
        self.bonds = nx.Graph()
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

    def _get_bond_fragments(atom1, atom2):
        '''
        Returns the pieces of a molecule that one would obtain by breaking the bond between two atoms. Will throw ValueError if the atoms are in a ring. 
        Useful for distance/angle/dihedral scans -- one fragment can be moved and the other held constant.   
        ''' 
        bond_order = self.bonds[atom1, atom2]['weight']
        if bond_order == 0:
            ValueError("No bond between atom {} and atom {}!".format(atom1,atom2))
        
        self.bonds.remove_edge(atom1, atom2)
        
        fragment1 = self.bonds[atom1]
        fragment2 = self.bonds[atom2]

        if atom1 in fragment2:
            ValueError("Atom {} and atom {} are in a ring or otherwise connected!".format(atom1,atom2)

        self.bonds.add_edge(atom1,atom2,weight=bond_order)

        return fragment1, fragment2


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
