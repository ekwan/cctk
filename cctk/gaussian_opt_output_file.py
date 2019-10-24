import sys
import re
import numpy as np

from cctk import GaussianOutputFile
from cctk.helper_functions import compute_distance_between, get_number

class GaussianOptOutputFile(GaussianOutputFile):
    """
    Class for Gaussian output files originating from a geometry optimization. 

    Attributes:
        atoms (list): list of atomic numbers 
        energies (list): list of energies for each cycle
        scf_iterations (list): number of iterations per cycle
        max_displacements (list): list of max displacement values for each cycle 
        rms_displacements (list): list of rms displacement values for each cycle 
        max_forces (list): list of max force values for each cycle 
        rms_forces (list): list of rms force values for each cycle 
        geometries (list): list of geometries for each cycle (so really a list of lists)
    """
    def __init__(self, filename):
        self.read_file(filename)   
        pass    

    def read_file(self, filename):
        """
        Reads a Gaussian optimization out file and populates the attributes accordingly. 

        Args:   
            filename (str): path to the out file
        """
        lines = super().read_file(filename)
        geometries, atom_list, energies, scf_iterations = super().read_geometries_and_energies(lines)

        self.atoms = list(map(get_number, atom_list))
        self.geometries = geometries
        self.energies = energies
        self.scf_iterations = scf_iterations
        
        self.rms_forces = self._find_parameter(lines, "RMS\s+Force")
        self.rms_displacements= self._find_parameter(lines, "RMS\s+Displacement")

    def _find_parameter(self, lines, parameter):
        """
        Helper method to search through the output file and find key forces and displacements. 

        Args:
            lines (list): list of lines from read_file()
            parameter (string): test to search for

        Returns:
            a list of all the extracted values
        """
        matches = []
        pattern = re.compile(parameter)
        for line in lines:
            if pattern.search(line):
                fields = re.split(' +', line)
                fields = list(filter(None, fields))
                if len(fields) != 5:
#                    print("error parsing >>> " + line)
#                    raise ValueError("unexpected number of fields on line")
                    continue
                matches.append(float(fields[2]))
        return matches

    def get_final_geometry(self):
        """
        Returns the last geometry from the out file.
        """
        return self.geometries[-1]
        
    def create_molecule():
        mol_geom = self.get_final_geometry()
        pass

    def create_molecule_from_geometry(self, number):
        pass 

    def print_geometric_parameters(self, parameter, atom1, atom2, atom3=None, atom4=None):
        """
        Computes and outputs geometric parameters (bond distances, angles, or dihedral angles) for every geometry. 

        Args:
            parameter (str): one of "angle", "distance", or "dihedral"
            atom1 (int): index (in self.atoms) of the atom in question
            atom2 (int): same, but for the second atom
            atom3 (int): same, but for the third atom (only required for parameter "angle" or "dihedral")
            atom4 (int): same, but for the fourth atom (only required for parameter "dihedral")
        
        Returns:
            a list of the specified parameter's values for each geometry
        """
        output = [None] * len(self.geometries)
        for index, geometry in enumerate(self.geometries):
            if parameter == "distance":
                output[index] = self.get_distance(geometry, atom1, atom2)
            elif parameter == "angle":
                pass
            elif parameter == "dihedral":
                pass
            else:
                ValueError("Invalid parameter {}!".format(parameter))
        
        return output            

    def get_distance(self, geometry, atom1, atom2):
        """
        Wrapper to compute distance between two atoms. 
        """
        v1 = np.array(geometry[atom1-1])
        v2 = np.array(geometry[atom2-1])
        return compute_distance_between(v1, v2)
