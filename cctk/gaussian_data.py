import numpy as np

from cctk import GaussianOutputFile, OutputFile
from cctk.helper_functions import (
    get_number,
    get_symbol,
    compute_distance_between,
    compute_angle_between,
    compute_dihedral_between,
)

import cctk.parse_gaussian as parse


class GaussianData(OutputFile):
    """
    Creates output file instances of the specific type through factory methods.
    """

    @classmethod
    def read_opt(cls, filename, return_lines=False):
        """
        Reads a Gaussian optimization out file and populates the attributes accordingly. 

        Args:   
            filename (str): path to the out file
            return_lines (Bool): whether the lines of the file should be returned

        Returns:
            GaussianOutputFile object
            (optional) the lines of the file
        """
        lines, header, success = cls._read_file(filename)
        (geometries, atom_list, energies, scf_iterations,) = parse.read_geometries_and_energies(lines)
        atoms = list(map(get_number, atom_list))
        bonds = parse.read_bonds(lines)

        f = GaussianOutputFile(atoms, geometries, bonds)

        f.energies = energies
        f.scf_iterations = scf_iterations
        f.header = header
        f.success = success

        f.rms_forces = parse.find_parameter(lines, "RMS\s+Force", expected_length=5, which_field=2)
        f.rms_displacements = parse.find_parameter(lines, "RMS\s+Displacement", expected_length=5, which_field=2)

        if return_lines:
            return f, lines
        else:
            return f

    @classmethod
    def read_opt_freq(cls, filename):
        """
        Reads the output of a Gaussian opt-freq job and returns the GaussianOutputFile object. 
        
        Args:   
            filename (str): path to the out file

        Returns:
            GaussianOutputFile object
        """
        f, lines = cls.read_opt(filename, return_lines=True)

        enthalpies = parse.find_parameter(lines, "thermal Enthalpies", expected_length=7, which_field=6)
        if len(enthalpies) == 1:
            f.enthalpy = enthalpies[0]
        elif len(enthalpies) > 1:
            raise ValueError("too many enthalpies found!")

        gibbs_vals = parse.find_parameter(lines, "thermal Free Energies", expected_length=8, which_field=7)
        if len(gibbs_vals) == 1:
            f.gibbs_free_energy = gibbs_vals[0]
        elif len(gibbs_vals) > 1:
            raise ValueError("too many gibbs free energies found!")

        frequencies = []
        try:
            frequencies += parse.find_parameter(lines, "Frequencies", expected_length=5, which_field=2)
            frequencies += parse.find_parameter(lines, "Frequencies", expected_length=5, which_field=3)
            frequencies += parse.find_parameter(lines, "Frequencies", expected_length=5, which_field=4)
            f.frequencies = sorted(frequencies)
        except:
            raise ValueError("error finding frequencies")

        return f

    @classmethod
    def _read_file(cls, filename):
        """
        Read a Gaussian output file. Automatically determines the theory line and if the job was successful. 
        
        Args:
            filename (str): path to the file

        Returns:
            lines (list): a list of all the lines
            header (str): the header string
            success (int): how many of the subjobs succeeded
        """
        lines = super().read_file(filename)
        header = parse.search_for_block(lines, "#p", "----")

        success = 0
        for line in lines:
            if line.strip().startswith("Normal termination"):
                success += 1

        return lines, header, success
