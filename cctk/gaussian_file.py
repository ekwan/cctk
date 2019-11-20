import sys
import re
import numpy as np

from enum import Enum

from cctk import File, Ensemble
from cctk.helper_functions import get_symbol, compute_distance_between, compute_angle_between, compute_dihedral_between, get_number

class job_type(Enum):
    """
    Class to contain allowed Gaussian job types. Not an exhaustive list, but should be fairly comprehensive.
    """
    SP = 0
    OPT = 1
    FREQ = 2
    IRC= 3
    NMR = 4
    POP = 5

class GaussianFile(File):
    """
    Class for Gaussian files. Composes ``Ensemble``.

    Attributes:
        molecules (Ensemble): ``Ensemble`` instance
        job_types (list): list of `job_type` instances
        header (str): optional, header of .gjf file
        footer (str): optional, footer of .gjf file
        success (int): number of successful terminations (should be 1 for an opt, 2 for opt and then freq, 1 for a single point energy, etc)
        theory (dict): contains information from header
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
        title (str): optional, title of .gjf file
    """

    def __init__(self, atoms, geometry, theory=None, header=None, footer=None, title="title", charge=0, multiplicity=1):
        """
        Create new GaussianInputFile object.
        """
        self.molecules = Ensemble(atoms, geometry, charge=charge, multiplicity=multiplicity)
        self.header = header
        self.footer = footer
        self.title = title

    def write_file(self, filename, memory=32, cores=16, chk_path=None):
        """
        Write a .gjf file, using object attributes.
        """

        if self.header:
            text = ""
            text += "%nprocshared={}GB\n".format(cores)
            text += "%mem={}GB\n".format(memory)

            if chk_path:
                text += "%chk={}\n".format(chk_path)

            text += self.header.rstrip()
            text += "\n"
            text += "\n"
            text += "{}\n".format(self.title)
            text += "\n"

            text += "{} {}\n".format(self.molecule.charge, self.molecule.multiplicity)
            for index, line in enumerate(self.molecule.geometry):
                text += "{:2d} {:.8f} {:.8f} {:.8f}\n".format(self.molecule.atoms[index], line[0], line[1], line[2])

            text += "\n"
            if self.footer:
                text += self.footer.rstrip()
                text += "\n"
                text += "\n"

            super().write_file(filename, text)
        else:
            raise ValueError("need header to write input file!")

    def num_imaginary(self):
        """
        Returns the number of imaginary frequencies.
        """
        return int(np.sum(np.array(self.frequencies) <= 0, axis=0))

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
