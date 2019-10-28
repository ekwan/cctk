import sys
import re
import numpy as np
import networkx as nx

from cctk.helper_functions import get_symbol, compute_distance_between, compute_unit_vector, get_covalent_radius

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
        
        if not all(isinstance(z, int) for z in atoms) or len(atoms) == 0:
            raise ValueError("invalid atom list")
        
        if not all(isinstance(x, tuple) for x in geometry) or len(geometry) == 0:
            raise ValueError("invalid geometry list")

        if not all(all(isinstance(y, float) for y in x) for x in geometry):
            raise TypeError("each element of self.geometry must be a 3-tuple")
       
        if not isinstance(theory, dict):
            raise TypeError("theory must be a dictionary!")
        
        if not isinstance(name, str):
            raise TypeError("name must be a string!")

        for atom in atoms: 
            try:
                get_symbol(atom)
            except ValueError:
                raise ValueError(f"unknwon atomic number {atom}")

        self.name = name
        self.atoms = atoms
        self.theory = theory
        self.geometry = geometry
        
        if bonds:
            pass
        else: 
            self.bonds = nx.Graph()
            self.bonds.add_nodes_from(range(1,len(atoms)+1))
            
    def assign_connectivity(self):            
        """
        Automatically recalculates bonds based on covalent radii. If two atoms are closer than the sum of their covalent raii + 0.5 Angstroms, then they are considered bonded. 
        """
            
        for i in range(1,len(self.atoms) + 1):
            for j in range(i+1, len(self.atoms) + 1):
                distance = self.get_distance(i, j)
                r_i = get_covalent_radius(self.get_atomic_number(i)) 
                r_j = get_covalent_radius(self.get_atomic_number(j)) 
                
                # 0.5 A distance is used by RasMol and Chime (documentation available online) and works well, empirically  
                if distance < (r_i + r_j + 0.5):
                    self.add_bond(i, j) 

    def add_bond(self, atom1, atom2, bond_order=1):
        """
        Adds a new bond to the bond graph, or updates the existing bond order. Will not throw an error if the bond already exists..

        Args:
            atom1 (int): the number of the first atom
            atom2 (int): the number of the second atom
            bond_order (int): bond order of bond between atom1 and atom2
        """
        self._check_atom_number(atom1)
        self._check_atom_number(atom2)

        if self.bonds.has_edge(atom1-1, atom2-1):
            if self.bonds[atom1][atom2]['weight'] != bond_order:
                self.bonds[atom1][atom2]['weight'] = bond_order
        else: 
            self.bonds.add_edge(atom1-1, atom2-1, weight=bond_order)

    def _check_atom_number(self, number):
        """
        Helper method which performs quick checks on the validity of a given atom number.
        """ 
        if not isinstance(number, int):
            raise TypeError("atom number must be integer")

        if number > len(self.atoms): 
            raise ValueError("atom numbers too large!")

        if atom1 <= 0:
            raise ValueError("atom number must be a positive integer!")
    
    def formula():
        pass

    def write_gaussian_input(header, footer=None):
        pass
    
    def append_to_gaussian_input(header, footer=None):
        '''
        Appends to GaussianInputFile object; to be appended in .gjf file using Link1 command
        '''
        pass

    def _get_bond_fragments(self, atom1, atom2, bond_order=1):
        '''
        Returns the pieces of a molecule that one would obtain by breaking the bond between two atoms. Will throw ValueError if the atoms are in a ring. 
        Useful for distance/angle/dihedral scans -- one fragment can be moved and the other held constant.   
        
        Args: 
            atom1 (int): the number of the first atom
            atom2 (int): the number of the second atom
            bond_order (int): bond order of bond between atom1 and atom2

        Returns:
            fragment1: the list of atoms in fragment 1 (containing atom1)
            fragment2: the list of atoms in fragment 2 (containing atom2)
        
        ''' 

        self._check_atom_number(atom1)
        self._check_atom_number(atom2)

        if (not isinstance(bond_order, int)) or (bond_order < 0):
            raise ValueError("invalid bond order!")

        # Users see indexing from one but system requires zero-indexing
        atom1 += -1
        atom2 += -1

        if self.bonds.has_edge(atom1, atom2):
            self.bonds.remove_edge(atom1, atom2)
            
            fragment1 = list(self.bonds[atom1].keys()) + [atom1]
            fragment2 = list(self.bonds[atom2].keys()) + [atom2]

            print(fragment1)
            print(fragment2)

            if atom1 in fragment2:
                raise ValueError(f"Atom {atom1+1} and atom {atom2+1} are in a ring or otherwise connected!")

            self.bonds.add_edge(atom1,atom2,weight=bond_order)
            return fragment1, fragment2
        else:
            raise ValueError(f"No bond between atom {atom1+1} and atom {atom2+1}!")
        
    def write_mol2_input(header, footer=None):
        pass

    def set_distance(self, atom1, atom2, distance, move='group'):
        """
        Adjusts the `atom1`&ndash;`atom2` bond length to be a fixed distance by moving atom2. 

        If `move` is set to "group", then all atoms bonded to `atom2` will also be moved. 
        
        If `move` is set to "atom", then only atom2 will be moved. 
        
        Args:
            atom1 (int): the number of the first atom
            atom2 (int): the number of the second atom
            distance (float): distance in Angstroms of the final bond
            move (str): determines how fragment moving is handled 
        """

        self._check_atom_number(atom1)
        self._check_atom_number(atom2)

        if (not isinstance(distance, float)) or (distance < 0):
            raise ValueError(f"invalid value {distance} for distance!")

        atoms_to_move = []
        if move == 'group':
            _, atoms_to_move = self._get_bond_fragments(atom1, atom2)
        elif move == 'atom':
            atoms_to_move = [atom2-1]
        else:
            raise ValueError(f"Invalid option {move} for parameter 'move'!")

        current_distance = self.get_distance(atom1,atom2)  
        
        v1 = self.get_geometry(atom1)
        v2 = self.get_geometry(atom2)
        vb = v2 - v1

        if np.linalg.norm(vb) - current_distance > 0.00001:
            raise ValueError(f"Error calculating bond distance!") 

        #### move all the atoms
        delta = distance - current_distance
        unitv = compute_unit_vector(vb) 
        for atom in atoms_to_move:
            self.geometry[atom] = self.geometry[atom] + (delta * unitv) 

        #### check everything worked okay...
        v1f = self.get_geometry(atom1)
        v2f = self.get_geometry(atom2)
        vbf = v2f - v1f
        
        if np.linalg.norm(vbf) - distance > 0.00001:
            raise ValueError(f"Error moving bonds -- operation failed!") 

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
        self_check_atom_number(atom)
        return self.atoms[atom-1]
    
    def get_geometry(self, atom):
        """
        Get the geometry vector for a given atom.
        """
        self_check_atom_number(atom)
        return np.array(self.geometry[atom-1])

    def get_distance(self, atom1, atom2):
        """
        Wrapper to compute distance between two atoms. 
        """
        self_check_atom_number(atom1)
        self_check_atom_number(atom2)

        v1 = np.array(self.geometry[atom1-1])
        v2 = np.array(self.geometry[atom2-1])
        return compute_distance_between(v1, v2)
