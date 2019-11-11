import sys
import re
import numpy as np

from cctk.helper_functions import align_matrices 

class Ensemble():
    """
    Class that represents a group of molecules.
    
    Attributes: 
        name (str): name, for identification
        molecules (list): list of `Molecule` objects
    """

    def __init__(self, name=None): 
        self.name = name
        self.molecules = []

    def align (self, align_to=1, atom_numbers=None):
        """ 
        Aligns every geometry to the specified geometry based on the atoms in `atom_numbers`. If `atom_numbers` is `None`, then a full alignment is performed. 
        
        Args: 
            align_to (int): which geometry to align to (1-indexed)
            atom_numbers (list): which atoms to align in each molecule (1-indexed; must be at least 3)
        """ 
        self._check_molecule_number(align_to)
        
        if atom_numbers and (len(atom_numbers) < 3):
            raise ValueError("not enough atoms for alignment - need 3 in 3D space!")

        try: 
            atom_numbers = np.array(atom_numbers)
            atom_numbers += -1
        except:
            raise ValueError("atom_numbers cannot be cast to numpy array... disappointing!")

        #### atom numbers is 0-indexed now
        #### move everything to the center!
        for molecule in self.molecules:
            centroid = self.geometry[atom_numbers].mean(axis=0)
            molecule.translate(-centroid)
        
        template = self.molecules[align_to-1].geometry[atom_numbers]
        
        #### perform alignment
        for molecule in self.molecules: 
            new_geometry = align_matrices(self.geometry[atom_numbers], self.geometry, template)
            self.geometry = new_geometry

            assert len(self.geometry) == len(self.atoms), "wrong number of geometry elements!"
        
    def eliminate_redundant(self, cutoff=0.05):
        pass    

    def lowest_energy(self, n=1):
        pass    

    def _check_molecule_number(self, number):
        """
        Helper method which performs quick checks on the validity of a given molecule number.
        """
        if not isinstance(number, int):
            raise TypeError("atom number must be integer")

        if number > len(self.moleculess):
            raise ValueError(f"atom number {number} too large!")

        if number <= 0:
            raise ValueError(f"atom number {number} invalid: must be a positive integer!")


