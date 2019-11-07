import sys
import re
import numpy as np

from cctk import OutputFile
from cctk.helper_functions import get_symbol, compute_distance_between, compute_angle_between, compute_dihedral_between, get_number

class GaussianOutputFile(OutputFile):
    '''
    Generic class for all Gaussian output files. 

    Attributes: 
        title (str): the title from the Gaussian file
        atoms (list): list of atomic numbers 
        geometries (list): list of geometries for each cycle (so really a list of lists, which may have only one element for single point jobs)
        bonds (list): list of atoms bonded to each other (so a list of 2-element lists)
        success (int): number of successful terminations (should be 1 for an opt, 2 for opt and then freq, 1 for a single point energy, etc)
        theory (dict): contains information from header
        header (str):  header from input file (e.g. "#p opt freq=noraman b3lyp/6-31g(d)").
        footer (str):  footer from input file, if any. genecp commands or opt=modredundant specifications go here. 
        energies (list): list of energies for each cycle
        scf_iterations (list): number of iterations per cycle
        max_displacements (list): list of max displacement values for each cycle 
        rms_displacements (list): list of rms displacement values for each cycle 
        max_forces (list): list of max force values for each cycle 
        rms_forces (list): list of rms force values for each cycle 
        gradients (list): list of gradient values for each cycle
        frequencies (list): list of frequencies
        gibbs_free_energy (float): gibbs free energy, from vibrational correction
        enthalpy (float): enthalpy, from vibrational correction
    '''
    def __init__(self, atoms, geometries, bonds=None):
        if (len(atoms) == 0) or (len(geometries) == 0): 
            raise ValueError("can't have a molecule without atoms or geometries!")

        for atom in atoms: 
            try:
                atom = int(atom)
                get_symbol(atom)
            except: 
                raise ValueError(f"atom number {atom} is not valid!")

        for geometry in geometries: 
            if (len(geometry) != len(atoms)):
                raise ValueError("length of atoms and geometry do not match!")        

        if bonds:
            for bond in bonds: 
                if len(bond) != 2:
                    raise ValueError("while 3-center bonding is possible, it's a no-go in cctk")
                if (not isinstance(bond[0], int)) or (not isinstance(bond[1], int)):
                    raise ValueError(f"atoms {bond[0]} and {bond[1]} must be integers!")
                if (bond[0] > len(atoms)) or (bond[0] <= 0):
                    raise ValueError(f"invalid atom {bond[0]}")
                if (bond[1] > len(atoms)) or (bond[1] <= 0):
                    raise ValueError(f"invalid atom {bond[1]}")
        
        self.atoms = atoms
        self.geometries = geometries
        self.bonds = bonds 

    def get_distance(self, geometry, atom1, atom2):
        """
        Wrapper to compute distance between two atoms. 
        """
        try:
            atom1 = int(atom1)
            atom2 = int(atom2)
        except:
            raise TypeError("atom1 and atom2 must be int!")
        
        if (atom1 <= len(geometry)) and (atom2 <= len(geometry)):
            v1 = np.array(geometry[atom1-1])
            v2 = np.array(geometry[atom2-1])
            return compute_distance_between(v1, v2)
        else:
            raise ValueError(f"atom numbers {atom1} or {atom2} too big!")

    def get_angle(self, geometry, atom1, atom2, atom3):
        """
        Wrapper to compute angle between two atoms. 
        """
        try:
            atom1 = int(atom1)
            atom2 = int(atom2)
            atom3 = int(atom3)
        except:
            raise TypeError("atom1, atom2, and atom3 must be int!")
        
        if ((atom1 <= len(geometry)) and (atom2 <= len(geometry))) and (atom3 <= len(geometry)):
            v1 = np.array(geometry[atom1-1])
            v2 = np.array(geometry[atom2-1])
            v3 = np.array(geometry[atom3-1])
            return compute_angle_between(v1, v2, v3)
        else:
            raise ValueError(f"atom numbers {atom1}, {atom2}, or {atom3} too big!")
   
    def get_dihedral(self, geometry, atom1, atom2, atom3, atom4):
        """
        Wrapper to compute dihedral angle between two atoms. 
        """
        try:
            atom1 = int(atom1)
            atom2 = int(atom2)
            atom3 = int(atom3)
            atom4 = int(atom4)
        except:
            raise TypeError("atom1, atom2, atom3, and atom4 must be int!")
        
        if ((atom1 <= len(geometry)) and (atom2 <= len(geometry))) and ((atom3 <= len(geometry)) and (atom4 <= len(geometry))):
            v1 = np.array(geometry[atom1-1])
            v2 = np.array(geometry[atom2-1])
            v3 = np.array(geometry[atom3-1])
            v4 = np.array(geometry[atom4-1])
            return compute_dihedral_between(v1, v2, v3, v4)
        else:
            raise ValueError(f"atom numbers {atom1}, {atom2}, {atom3}, or {atom4} too big!")

    def num_imaginary(self):
        """ 
        Returns the number of imaginary frequencies.
        """
        return int(np.sum(np.array(self.frequencies) <= 0, axis=0))
    
    def get_final_geometry(self):
        """
        Returns the last geometry from the out file.
        """
        return self.geometries[-1]
        
    def print_geometric_parameters(self, parameter, atom1, atom2, atom3=None, atom4=None):
        """
        Computes and outputs geometric parameters (bond distances, angles, or dihedral angles) for every geometry. 

        Args:
            parameter (str): one of ``angle``, ``distance``, or ``dihedral``
            atom1 (int): number of the atom in question 
            atom2 (int): same, but for the second atom
            atom3 (int): same, but for the third atom (only required for parameter ``angle`` or ``dihedral``)
            atom4 (int): same, but for the fourth atom (only required for parameter ``dihedral``)
        
        Returns:
            a list of the specified parameter's values for each geometry
        """
        output = [None] * len(self.geometries)
        for index, geometry in enumerate(self.geometries):
            if parameter == "distance":
                output[index] = self.get_distance(geometry, atom1, atom2)
            elif parameter == "angle":
                if atom3 == None: 
                    raise ValueError("need atom3 to calculate angle!")
                output[index] = self.get_angle(geometry, atom1, atom2, atom3)
            elif parameter == "dihedral":
                if (atom3 == None) or (atom4 == None):
                    raise ValueError("need atom3 and atom4 to calculate dihedral!")
                output[index] = self.get_dihedral(geometry, atom1, atom2, atom3, atom4)
            else:
                ValueError("Invalid parameter {}!".format(parameter))
        
        return output            

