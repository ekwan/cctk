import sys
import re
import numpy as np
import networkx as nx

from cctk.helper_functions import get_symbol, compute_rotation_matrix, compute_distance_between, compute_angle_between, compute_dihedral_between, compute_unit_vector, get_covalent_radius

class Molecule():
    """
    Class that represents a single molecule, abstractly.

    Outward-facing methods are indexed from one (i.e. ``self.get_vector(1)`` returns the first element of self.geometry). Internal methods (prefixed with ``_``) are generally zero-indexed, and list attributes (e.g. ``self.geometry``, ``self.atoms``) are zero-indexed inherently. 

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
       
        if len(geometry) == 0:
            raise ValueError("invalid geometry list")

        try:
            geometry = np.array(geometry)
            geometry = geometry.astype(float)
        except:
            raise TypeError("geometry cannot be cast to ndarray of floats!")  

        if not all(all(isinstance(y, float) for y in x) for x in geometry):
            raise TypeError("each element of self.geometry must be a 3-tuple")
       
        if theory and not isinstance(theory, dict):
            raise TypeError("theory must be a dictionary!")
        
        if name and not isinstance(name, str):
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
        
        self.bonds = nx.Graph()
        self.bonds.add_nodes_from(range(1,len(atoms)+1))
        if bonds:
            for bond in bonds:
                self._check_atom_number(bond[0])
                self._check_atom_number(bond[1])
                self.add_bond(bond[0], bond[1]) 
            
    def assign_connectivity(self):            
        """
        Automatically recalculates bonds based on covalent radii. If two atoms are closer than the sum of their covalent radii + 0.5 Angstroms, then they are considered bonded. 
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
            raise ValueError(f"atom number {number} too large!")

        if number <= 0:
            raise ValueError(f"atom number {number} invalid: must be a positive integer!")
    
    def formula(self, string=False):
        """
        Returns the atomic formula. 
        
        If ``string`` is ``False``, then returns a ``dictionary`` with keys elemental symbols and values the number of occurrences. 
        
        For instance, ``water.formula()`` would return ``{'O': 1, 'H': 2}``. 
        
        If ``string`` is ``True``, then returns a stringified version of the formula according to standard rules. 

        For instance, ``water.formula()`` would return ``H2O``. 
        
        Args:
            string (Bool): if the method should return a string or a dictionary

        Returns:
            a dictionary or string representing the molecule's formula
        """

        formula_dict = {}
        for atom in self.atoms:
            symbol = get_symbol(atom)
            if symbol in formula_dict:
                formula_dict[symbol] += 1
            else:
                formula_dict[symbol] = 1
        if string == False:
            return formula_dict
        else:
            formula = '' 
            elements = list(formula_dict.keys())
                    
            #### H and C always come first  
            if 'H' in elements: 
                elements.remove('H')
                formula += f"H{formula_dict['H']}"

            if 'C' in elements: 
                elements.remove('C')
                formula += f"C{formula_dict['C']}"

            for element in sorted(elements):
                formula += f"{element}{formula_dict[element]}"
            
            return formula 

    def _get_bond_fragments(self, atom1, atom2, bond_order=1):
        '''
        Returns the pieces of a molecule that one would obtain by breaking the bond between two atoms. Will throw ``ValueError`` if the atoms are in a ring. 
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
           
            fragments = nx.connected_components(self.bonds)
            fragment1 = []
            fragment2 = []

            for fragment in fragments:
                if atom1 in fragment:
                    if atom2 in fragment:
                        raise ValueError(f"Atom {atom1+1} and atom {atom2+1} are in a ring or otherwise connected!")
                    else:
                        fragment1 = fragment
                if atom2 in fragment:
                    fragment2 = fragment

            self.bonds.add_edge(atom1,atom2,weight=bond_order)
            return fragment1, fragment2
        else:
            raise ValueError(f"No bond between atom {atom1+1} and atom {atom2+1}!")
       
    def _get_fragment_containing(self, atom):
        """ 
        Get the fragment containing the atom with number ``atom``.
        
        Args:
            atom (int): the number of the atom (indexed from one)
        
        Returns:
            a (zero-indexed) list of all the atoms in the fragment 
        """

        self._check_atom_number(atom)
        fragments = nx.connected_components(self.bonds)
        
        #### fragment is zero-indexed
        for fragment in fragments:
            if atom-1 in fragment: 
                return fragment
                break

    def set_distance(self, atom1, atom2, distance, move='group'):
        """
        Adjusts the ``atom1`` -- ``atom2`` bond length to be a fixed distance by moving atom2. 

        If ``move`` is set to "group", then all atoms bonded to ``atom2`` will also be moved. 
        
        If ``move`` is set to "atom", then only atom2 will be moved. 
        
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
            if self.get_bond_order(atom2, atom3):
                _, atoms_to_move = self._get_bond_fragments(atom2, atom3)
            else:
                atoms_to_move = self._get_fragment_containing(atom3)
        elif move == 'atom':
            atoms_to_move = [atom2-1]
        else:
            raise ValueError(f"Invalid option {move} for parameter 'move'!")

        current_distance = self.get_distance(atom1,atom2)  
        
        v1 = self.get_vector(atom1)
        v2 = self.get_vector(atom2)
        vb = v2 - v1

        if np.linalg.norm(vb) - current_distance > 0.00001:
            raise ValueError(f"Error calculating bond distance!") 

        #### move all the atoms
        delta = distance - current_distance
        unitv = compute_unit_vector(vb) 
        for atom in atoms_to_move:
            self.geometry[atom] = self.geometry[atom] + (delta * unitv) 

        #### check everything worked okay...
        v1f = self.get_vector(atom1)
        v2f = self.get_vector(atom2)
        vbf = v2f - v1f
        
        if np.linalg.norm(vbf) - distance > 0.001:
            raise ValueError(f"Error moving bonds -- operation failed!") 

    def set_angle(self, atom1, atom2, atom3, angle, move='group'):
        """
        Adjusts the ``atom1`` -- ``atom2`` -- ``atom3`` bond angle to be a fixed value by moving ``atom3``. 

        If `move` is set to "group", then all atoms bonded to ``atom3`` will also be moved. 
        
        If `move` is set to "atom", then only ``atom3`` will be moved. 
        
        Args:
            atom1 (int): the number of the first atom
            atom2 (int): the number of the second atom
            atom3 (int): the number of the third atom
            angle (float): final value in degrees of the ``atom1`` -- ``atom2`` -- ``atom3`` angle
            move (str): determines how fragment moving is handled 
        """

        self._check_atom_number(atom1)
        self._check_atom_number(atom2)
        self._check_atom_number(atom3)

        if self.get_distance(atom1, atom2) < 0.01: 
            raise ValueError(f"atom {atom1} and atom {atom2} are too close!")

        if self.get_distance(atom2, atom3) < 0.01: 
            raise ValueError(f"atom {atom2} and atom {atom3} are too close!")
        
        if self.get_distance(atom1, atom3) < 0.01: 
            raise ValueError(f"atom {atom1} and atom {atom3} are too close!")

        try:
            angle = float(angle)
        except: 
            raise TypeError(f"angle {angle} cannot be converted to float!")

        if (not isinstance(angle, float)) or ((angle < 0) or (angle > 360)):
            raise ValueError(f"invalid value {angle} for angle!")

        atoms_to_move = []
        if move == 'group':
            if self.get_bond_order(atom2, atom3):
                _, atoms_to_move = self._get_bond_fragments(atom2, atom3)
            elif self.are_connected(atom2, atom3):
                raise ValueError(f"atom {atom2} and atom {atom3} are connected but not bonded -- cannot adjust angle! try manually removing one or more bonds.")
            else:
                atoms_to_move = self._get_fragment_containing(atom3)
        elif move == 'atom':
            atoms_to_move = [atom3-1]
        else:
            raise ValueError(f"Invalid option {move} for parameter 'move'!")

        #### atoms to move is zero-indexed
        if atom1-1 in atoms_to_move:
            raise ValueError(f"atom {atom1} and atom {atom3} are connected in multiple ways -- cannot adjust angle! try manually removing one or more bonds.")

        current_angle = self.get_angle(atom1, atom2, atom3)
        delta = angle - current_angle

        if np.abs(delta) < 0.001: 
            return

        #### now the real work begins...

        #### move everything to place atom2 at the origin
        v2 = self.get_vector(atom2)
        self.translate_molecule(-v2)        

        v1 = self.get_vector(atom1)
        v3 = self.get_vector(atom3)
        
        #### perform the actual rotation
        rot_axis = np.cross(v1,v3)
        rot_matrix = compute_rotation_matrix(rot_axis, delta) 

        for atom in atoms_to_move:
            #### have to add one because atoms_to_move is zero indexed while get_vector is one indexed
            self.geometry[atom] = list(np.dot(rot_matrix, self.get_vector(atom+1)))
        
        #### and move it back!
        self.translate_molecule(v2) 
        
        final_angle = self.get_angle(atom1, atom2, atom3)
        if np.abs(final_angle - angle) > 0.001:
            raise ValueError(f"Error rotating atoms -- expected angle {angle}, got {final_angle}  -- operation failed!")

    def set_dihedral(self, atom1, atom2, atom3, atom4, dihedral, move='group34'):
        """
        Adjusts the ``atom1`` -- ``atom2`` -- ``atom3`` -- ``atom4`` dihedral angle to be a fixed value by moving atom 4.. 

        If ``move`` is set to "atom", then only ``atom4`` will be moved. 
        
        If ``move`` is set to "group4", then all atoms bonded to ``atom4`` will also be moved. 
       
        If ``move`` is set to "group34", then all atoms bonded to ``atom3`` and ``atom4`` will also be moved. 
        
        Args:
            atom1 (int): the number of the first atom
            atom2 (int): the number of the second atom
            atom3 (int): the number of the third atom
            atom4 (int): the number of the fourth atom
            dihedral (float): final value in degrees of the ``atom1`` -- ``atom2`` -- ``atom3`` -- ``atom4`` angle
            move (str): determines how fragment moving is handled 
        """

        self._check_atom_number(atom1)
        self._check_atom_number(atom2)
        self._check_atom_number(atom3)
        self._check_atom_number(atom4)

        for x in [atom1, atom2, atom3, atom4]:
            for y in [atom1, atom2, atom3, atom4]:
                if x == y:
                    continue
                else: 
                    if self.get_distance(x, y) < 0.01:
                        raise ValueError(f"atom {x} and atom {y} are too close!")
        
        try:
            dihedral = float(dihedral)
        except: 
            raise TypeError(f"dihedral angle {dihedral} cannot be converted to float!")

        if (not isinstance(dihedral, float)) or ((dihedral < 0) or (dihedral > 360)):
            raise ValueError(f"invalid value {dihedral} for dihedral angle!")

        atoms_to_move = []
        if move == 'group34':
            #### add atom3's fragment to atom4
            if self.get_bond_order(atom2, atom3):
                _, atoms_to_move = self._get_bond_fragments(atom2, atom3)
            elif self.are_connected(atom2, atom3):
                raise ValueError(f"atom {atom2} and atom {atom3} are connected but not bonded -- cannot adjust dihedral angle! try manually removing one or more bonds.")
            else:
                atoms_to_move = self._get_fragment_containing(atom3)
            
            #### and make sure atom4 is in there too!
            if atom4-1 not in atoms_to_move:
                atoms_to_move += self._get_fragment_containing(atom4)
        elif move == 'group4':
            if self.get_bond_order(atom3, atom4):
                _, atoms_to_move = self._get_bond_fragments(atom3, atom4)
            elif self.are_connected(atom3, atom4):
                raise ValueError(f"atom {atom3} and atom {atom4} are connected but not bonded -- cannot adjust dihedral angle! try manually removing one or more bonds.")
            else:
                atoms_to_move = self._get_fragment_containing(atom4)
        elif move == 'atom':
            atoms_to_move = [atom4-1]
        else:
            raise ValueError(f"Invalid option {move} for parameter 'move'!")

        #### atoms to move is zero-indexed
        if atom1-1 in atoms_to_move:
            raise ValueError(f"atom {atom1} and atom {atom4} are connected in multiple ways -- cannot adjust dihedral angle! try manually removing one or more bonds.")
        
        if atom2-1 in atoms_to_move:
            raise ValueError(f"atom {atom2} and atom {atom4} are connected in multiple ways -- cannot adjust dihedral angle! try manually removing one or more bonds.")
        
        if atom4-1 not in atoms_to_move:
            raise ValueError(f"atom {atom4} is not going to be moved... this operation is doomed to fail!")
        
        current_dihedral = self.get_dihedral(atom1, atom2, atom3, atom4)
        delta = dihedral - current_dihedral

        if np.abs(delta) < 0.001: 
            return

        #### now the real work begins...

        #### move everything to place atom2 at the origin
        v3 = self.get_vector(atom3)
        self.translate_molecule(-v3)        

        v1 = self.get_vector(atom1)
        v2 = self.get_vector(atom2)
        v4 = self.get_vector(atom4)
        
        #### perform the actual rotation
        #### rot_axis is formally v3 - v2, but v3 is the zero vector after the translation... 
        rot_axis = -v2 
        rot_matrix = compute_rotation_matrix(rot_axis, delta) 

        for atom in atoms_to_move:
            #### have to add one because atoms_to_move is zero indexed while get_vector is one indexed
            self.geometry[atom] = list(np.dot(rot_matrix, self.get_vector(atom+1)))
        
        #### and move it back!
        self.translate_molecule(v3) 
        
        final_dihedral = self.get_dihedral(atom1, atom2, atom3, atom4)
        if np.abs(final_dihedral - dihedral) > 0.001:
            raise ValueError(f"Error rotating atoms -- expected dihedral angle {dihedral}, got {final_dihedral}  -- operation failed!")

    def translate_molecule(self, vector):
        for atom in range(0,len(self.atoms)): 
            self.geometry[atom] = list(self.geometry[atom] + vector)

    def rotate_molecule(self, axis, degrees):
        pass

    def calculate_mass_spectrum():
        pass
    
    def add_atom_at_centroid(self, symbol, atom_numbers):
        pass

    def get_atomic_number(self, atom):
        """ 
        Get the atomic number for a given atom.
        """
        self._check_atom_number(atom)
        return self.atoms[atom-1]
    
    def get_vector(self, atom):
        """
        Get the geometry vector for a given atom.
        """
        self._check_atom_number(atom)
        return np.array(self.geometry[atom-1])

    def get_distance(self, atom1, atom2):
        """
        Wrapper to compute distance between two atoms. 
        """
        self._check_atom_number(atom1)
        self._check_atom_number(atom2)

        return compute_distance_between(self.get_vector(atom1), self.get_vector(atom2))

    def get_angle(self, atom1, atom2, atom3):
        """
        Wrapper to compute angle between three atoms. 
        """
        self._check_atom_number(atom1)
        self._check_atom_number(atom2)
        self._check_atom_number(atom3)
       
        v1 = self.get_vector(atom1)
        v2 = self.get_vector(atom2)
        v3 = self.get_vector(atom3)
        
        return compute_angle_between(v1-v2, v3-v2)
     
    def get_dihedral(self, atom1, atom2, atom3, atom4):
        """
        Wrapper to compute angle between three atoms. 
        """
        self._check_atom_number(atom1)
        self._check_atom_number(atom2)
        self._check_atom_number(atom3)
        self._check_atom_number(atom4)
       
        v1 = self.get_vector(atom1)
        v2 = self.get_vector(atom2)
        v3 = self.get_vector(atom3)
        v4 = self.get_vector(atom4)
        
        return compute_dihedral_between(v1, v2, v3, v4)


    def get_bond_order(self, atom1, atom2):
        """
        Wrapper to get bond order between two atoms.
        """
        self._check_atom_number(atom1)
        self._check_atom_number(atom2)
        
        if self.bonds.has_edge(atom1-1, atom2-1):
            return self.bonds[atom1-1][atom2-1]['weight']
        else:
            return 0

    def are_connected(self, atom1, atom2):
        """
        Wrapper to tell if two atoms are connected. 
        """
        self._check_atom_number(atom1)
        self._check_atom_number(atom2)
        
        if atom1-1 in self._get_fragment_containing(atom2):
            return True
        else:
            return False
