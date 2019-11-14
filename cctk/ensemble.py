import sys
import re
import numpy as np
import copy

from cctk import Molecule
from cctk.helper_functions import align_matrices, compute_RMSD

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
        
    def eliminate_redundant(self, cutoff=0.5, heavy_only=True, atom_numbers=None):
        """
        Returns non-redundant conformations. When redundancies are found, only the first geometry is kept.
        This will change the numbering of all the ensembles!
        
        Args:
            cutoff (float): molecules with less than this value for RMSD will be considered redundant and eliminated. 
            heavy_only (Bool): if `True`, then only heavy atoms are considered for the RMSD calculation
            atom_numbers (list): 1-indexed list of atoms to consider for RMSD calculation - if present, overrides `heavy_only`
        """
        if atom_numbers:
            atom_numbers = [n-1 for n in atom_numbers]
        else:
            if heavy_only: 
                atom_numbers = self.molecules[0].get_heavy_atoms()        
            else:
                atom_numbers = list(range(len(self.molecules[0])))
        
        for m in self.molecules:
            for n in atom_numbers:
                try:
                    #### atom_numbers is 0-indexed
                    m._check_atom_number(n+1)
                except:
                    raise ValueError(f"molecule in ensemble does not have atom {n}!")
       
        #### align all molecules 
        self.align(atoms=atom_numbers)

        to_delete = [False] * len(self.molecules)

        for i in range(len(self.molecules)):
            if to_delete[i]:
                continue
            for j in range(i+1, len(self.molecules)):
                if to_delete[j]:
                    continue
                
                geometry1 = self.molecules[i].geometry[atom_numbers]
                geometry2 = self.molecules[j].geometry[atom_numbers]
              
                rmsd = compute_RMSD(geometry1, geometry2) 
                if rmsd < cutoff:
                    to_delete[j] = True       

        #### you have to delete in reverse order or you'll throw off the subsequent indices 
        for i in sorted(range(len(self.molecules)), reverse=True):
            if to_delete[i]:
                del self.molecules[i]

    def _check_molecule_number(self, number):
        """
        Helper method which performs quick checks on the validity of a given molecule number.
        """
        try:
            number = int(number)
        except:
            raise TypeError(f"atom number {number} must be integer")

        if number > len(self.molecules):
            raise ValueError(f"atom number {number} too large!")

        if number <= 0:
            raise ValueError(f"atom number {number} invalid: must be a positive integer!")


