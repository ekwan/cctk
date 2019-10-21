import sys
import re
import numpy as np

from cctk import GaussianOutputFile

class GaussianOptOutputFile(GaussianOutputFile):
    '''
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
    def get_lowest_energy_geometry():
        pass
        
    def create_molecule():
        mol_geom = self.get_lowest_energy_geometry()
        pass

    def create_molecule_from_geometry(number):
        pass 

    def print_energies():
        pass

    def print_geometric_parameters():
        for geometry in self.geometries():
            pass
