import sys
import re
import numpy as np

from cctk import InputFile
from cctk.helper_functions import get_symbol, compute_distance_between, compute_angle_between, compute_dihedral_between, get_number


class XYZFile(InputFile):
    """
    Generic class for all xyz files. 

    Attributes: 
        title (str): the title from the file
        atoms (list): list of atomic numbers 
        geometries (list): list of 3-element lists denoting the geometry of each point 
    """

    def __init__(self, atoms, geometry, title=None):
        if (len(atoms) == 0) or (len(geometry) == 0):
            raise ValueError("can't have a molecule without atoms or geometry!")

        for atom in atoms:
            try:
                atom = int(atom)
                get_symbol(atom)
            except:
                raise ValueError(f"atom number {atom} is not valid!")

        if len(geometry) != len(atoms):
            raise ValueError("length of atoms and geometry do not match!")

        if title and (isinstance(title, str)):
            self.title = title

        self.atoms = atoms
        self.geometry = geometry

    def get_distance(self, atom1, atom2):
        """
        Wrapper to compute distance between two atoms. 
        """
        try:
            atom1 = int(atom1)
            atom2 = int(atom2)
        except:
            raise TypeError("atom1 and atom2 must be int!")

        if (atom1 <= len(self.geometry)) and (atom2 <= len(self.geometry)):
            v1 = np.array(self.geometry[atom1 - 1])
            v2 = np.array(self.geometry[atom2 - 1])
            return compute_distance_between(v1, v2)
        else:
            raise ValueError(f"atom numbers {atom1} or {atom2} too big!")

    def get_angle(self, atom1, atom2, atom3):
        """
        Wrapper to compute angle between two atoms. 
        """
        try:
            atom1 = int(atom1)
            atom2 = int(atom2)
            atom3 = int(atom3)
        except:
            raise TypeError("atom1, atom2, and atom3 must be int!")

        if ((atom1 <= len(self.geometry)) and (atom2 <= len(self.geometry))) and (atom3 <= len(self.geometry)):
            v1 = np.array(self.geometry[atom1 - 1])
            v2 = np.array(self.geometry[atom2 - 1])
            v3 = np.array(self.geometry[atom3 - 1])
            return compute_angle_between(v1, v2, v3)
        else:
            raise ValueError(f"atom numbers {atom1}, {atom2}, or {atom3} too big!")

    def get_dihedral(self, atom1, atom2, atom3, atom4):
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

        if ((atom1 <= len(self.geometry)) and (atom2 <= len(self.geometry))) and ((atom3 <= len(self.geometry)) and (atom4 <= len(self.geometry))):
            v1 = np.array(self.geometry[atom1 - 1])
            v2 = np.array(self.geometry[atom2 - 1])
            v3 = np.array(self.geometry[atom3 - 1])
            v4 = np.array(self.geometry[atom4 - 1])
            return compute_dihedral_between(v1, v2, v3, v4)
        else:
            raise ValueError(f"atom numbers {atom1}, {atom2}, {atom3}, or {atom4} too big!")

    def write_file(self, filename):
        """ 
        Write a .xyz file, using object attributes. "title" will be used as the title if none is defined.  
        
        Args:
            filename (str): path to the new file
        """
        text = f"{len(self.atoms)}\n"
        if self.title:
            text += f"{self.title}\n"
        else:
            text += "title\n"
        for index, line in enumerate(self.geometry):
            text += "{:2s} {:.8f} {:.8f} {:.8f}\n".format(get_symbol(self.atoms[index]), line[0], line[1], line[2])

        super().write_file(filename, text)
