import sys
import re
import numpy as np
import networkx as nx

from cctk.helper_functions import compute_distance_between, get_covalent_radius

class Molecule():
    """
    Class that represents a single molecule, abstractly.

    Attributes:
        name (str): for identification, optional
        theory (dict): containing information about theory, to keep track for basis set scans, optional
        atoms (list): list of atomic symbols 
        geometry (list):  list of 3-tuples of xyz coordinates
        bonds (nx.Graph): Graph object containing connectivity information
        energy (float): holds the calculated energy of the molecule, optional
    """
    
    def __init__(self, atoms, geometry, name=None, theory=None, bonds=None):
        """
        Create new Molecule object, and assign connectivity if needed. 
        """
        if len(atoms) != len(geometry):
            raise ValueError("length of geometry and atoms does not match!")
        
        self.name = name
        self.atoms = atoms
        self.theory = theory
        self.geometry = geometry
        if bonds:
            pass
        else: 
            self.bonds = nx.Graph()
            
    def assign_connectivity(self):            
        for i in range(1,len(self.atoms) + 1):
            for j in range(i, len(self.atoms) + 1):
                distance = self.get_distance(i, j)
                r_i = get_covalent_radius(self.get_atomic_number(i)) 
                r_j = get_covalent_radius(self.get_atomic_number(j)) 
                
                if distance < (r_i + r_j):
                    self.add_bond(i, j) 

    def add_bond(self, atom1, atom2):
        """
        Adds a new bond to the bonding graph.

        Args:
            atom1 (int): the number of the first atom
            atom2 (int): the number of the second atom
        """
        self.bonds.add_edge(atom1-1, atom2-1)

    def formula():
        pass

    def write_gaussian_input(header, footer=None):
        pass
    
    def append_to_gaussian_input(header, footer=None):
        '''
        Appends to GaussianInputFile object; to be appended in .gjf file using Link1 command
        '''
        pass

    def _get_bond_fragments(self, atom1, atom2):
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
            ValueError("Atom {} and atom {} are in a ring or otherwise connected!".format(atom1,atom2))

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

    def add_atom_at_centroid(self, symbol, atom_numbers):
        pass

    def get_atomic_number(self, atom):
        """ 
        Get the atomic number for a given atom.
        """
        return self.atoms[atom-1]
    
    def get_distance(self, atom1, atom2):
        """
        Wrapper to compute distance between two atoms. 
        """
        v1 = np.array(self.geometry[atom1-1])
        v2 = np.array(self.geometry[atom2-1])
        return compute_distance_between(v1, v2)
