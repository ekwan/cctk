import sys
import re
import numpy as np
import copy

from cctk import Molecule
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

    def add_molecule(self, molecule):
        """
        Adds a molecule to the ensemble. `copy.deepcopy` is used so that an independent copy of the molecule is saved. 

        Args:
            molecule (Molecule): the molecule to be added
        """
        if not isinstance(molecule, Molecule):
            raise TypeError("molecule is not a Molecule - so it can't be added!")

        self.molecules.append(copy.deepcopy(molecule))

    def align (self, align_to=1, atoms=None):
        """ 
        Aligns every geometry to the specified geometry based on the atoms in `atom_numbers`. If `atom_numbers` is `None`, then a full alignment is performed. 
        
        Args: 
            align_to (int): which geometry to align to (1-indexed)
            atoms (list): which atoms to align in each molecule (1-indexed; must be at least 3)
        """ 
        self._check_molecule_number(align_to)
        
        if atoms and (len(atoms) < 3):
            raise ValueError("not enough atoms for alignment - need 3 in 3D space!")

        try: 
            atoms = np.array(atoms)
            atoms += -1
        except:
            raise ValueError("atom_numbers cannot be cast to numpy array... disappointing!")

        #### atom numbers is 0-indexed now
        #### move everything to the center!
        for molecule in self.molecules:
            centroid = molecule.geometry[atoms].mean(axis=0)
            molecule.translate_molecule(-centroid)
        
        template = self.molecules[align_to-1].geometry[atoms]
        
        #### perform alignment
        for molecule in self.molecules: 
            new_geometry = align_matrices(molecule.geometry[atoms], molecule.geometry, template)
            self.geometry = new_geometry

            assert len(molecule.geometry) == len(molecule.atoms), "wrong number of geometry elements!"
        
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

        if number > len(self.molecules):
            raise ValueError(f"atom number {number} too large!")

        if number <= 0:
            raise ValueError(f"atom number {number} invalid: must be a positive integer!")


