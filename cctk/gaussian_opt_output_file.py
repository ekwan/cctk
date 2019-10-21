import sys
import re
import numpy as np

from cctk import GaussianOutputFile
from cctk.helper_functions import compute_distance_between

class GaussianOptOutputFile(GaussianOutputFile):
    '''
    atoms = list of atoms 
    energies = list of energies for each cycle
    scf_iterations = number of iterations per cycle
    rms_forces
    max_forces
    rms_displacements
    max_displacements
    geometries = list of geometrries for each cycle
    done = Boolean
    '''
    def __init__(self, filename):
        self.read_file(filename)   
        pass    

    def read_file(self, filename):
        lines = super().read_file(filename)
        geometries, atom_list, energies, scf_iterations = super().read_geometries_and_energies(lines)

        self.atoms = atom_list
        self.geometries = geometries
        self.energies = energies
        self.scf_iterations = scf_iterations
        
        self.rms_forces = self._find_parameter(lines, "RMS\s+Force")
        self.rms_displacements= self._find_parameter(lines, "RMS\s+Displacement")

    def _find_parameter(self, lines, parameter):
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
        return geometries[-1]
        
    def create_molecule():
        mol_geom = self.get_final_geometry()
        pass

    def create_molecule_from_geometry(self, number):
        pass 

    def print_geometric_parameters(self, parameter, atom1, atom2, atom3=None, atom4=None):
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
        v1 = np.array(geometry[atom1-1])
        v2 = np.array(geometry[atom2-1])
        return compute_distance_between(v1, v2)
