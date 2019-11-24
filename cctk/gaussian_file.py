import sys
import re
import numpy as np

from enum import Enum

from cctk import File, Ensemble
from cctk.helper_functions import get_symbol, compute_distance_between, compute_angle_between, compute_dihedral_between, get_number

import cctk.parse_gaussian as parse

class JobType(Enum):
    """
    Class to contain allowed Gaussian job types. Not an exhaustive list, but should be fairly comprehensive.

    The value should be the Gaussian keyword, to permit automatic assignment.
    """
    SP = "sp"
    OPT = "opt"
    FREQ = "freq"
    IRC= "irc"
    NMR = "nmr"
    POP = "pop"

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

    def __init__(self, atomic_numbers, geometries, bonds=None, job_types=None, theory=None, header=None, footer=None, title="title", charge=0, multiplicity=1):
        """
        Create new GaussianInputFile object.

		Args:
            atomic_numbers (list): list of atomic numbers
            geometries (list): list of lists of 3-tuples of xyz coordinates
            bonds (nx.Graph): Graph object containing connectivity information (1-indexed)
            charge (int): the charge of the molecule
            multiplicity (int): the spin state of the molecule (1 corresponds to singlet, 2 to doublet, 3 to triplet, etc. -- so a multiplicity of 1 is equivalent to S=0)
            job_types (list): list of `job_type` instances
            header (str): optional, header of .gjf file
            footer (str): optional, footer of .gjf file
            theory (dict): optional, contains information from header
            title (str): optional, title of .gjf file
		"""

        if header and not isinstance(header, str):
            raise TypeError("header needs to be a string")

        if footer and not isinstance(footer, str):
            raise TypeError("footer needs to be a string")

        if title and not isinstance(title, str):
            raise TypeError("title needs to be a string")

        if not all(isinstance(job, JobType) for job in job_types):
            raise TypeError(f"invalid job type {job}")

        self.molecules = Ensemble(atomic_numbers=atomic_numbers, geometries=geometries, bonds=bonds, charge=charge, multiplicity=multiplicity)
        self.header = header
        self.footer = footer
        self.title = title
        self.job_types = job_types

    def write_file(self, filename, memory=32, cores=16, chk_path=None, molecule=None, header=None, footer=None):
        """
        Write a .gjf file, using object attributes. If no header is specified, the object's header/footer will be used.

        Args:
            filename (str): path to the new file
            memory (int): how many GB of memory to request
            cores (int): how many CPU cores to request
            chk_path (str): path to checkpoint file, if desired
            molecule (int): which molecule to use -- passed to ``self.get_molecule()``. Default is -1 (e.g. the last molecule), but positive integers will select from self.molecules (1-indexed).
            header (str): header for new file
            footer (str): footer for new file
        """
        mol = self.get_molecule(molecule)

        if header is None:
            header = self.header
            footer = self.footer

        #### generate the text
        text = f"%nprocshared={int(cores)}GB\n"
        text += f"%mem={int(memory)}GB\n"

        if chk_path:
            text += f"%chk={chk_path}\n"

        text += f"{header.strip()}\n\n{self.title}\n\n"

        text += f"{int(mol.charge)} {int(mol.multiplicity)}\n"
        for index, line in enumerate(mol.geometry):
            text += f"{mol.atomic_numbers[index]:2d} {line[0]:.8f} {line[1]:.8f} {line[2]:.8f}\n"

        text += "\n"
        if footer is not None:
            text += f"{footer.strip()}\n\n"

        #### write the file
        super().write_file(filename, text)

    def num_imaginary(self):
        """
        Returns the number of imaginary frequencies.
        """
        if JobType.FREQ in self.job_types:
            return int(np.sum(np.array(self.frequencies) <= 0, axis=0))
        else:
            raise TypeError("not a frequency job! can't get # imaginary frequencies!")

    @classmethod
    def read_file(cls, filename, job_types=[], return_lines=False):
        """
        Reads a Gaussian optimization out file and populates the attributes accordingly.

        Args:
            filename (str): path to the out file
            return_lines (Bool): whether the lines of the file should be returned
            job_types (list): list of JobTypes - if None, will be automatically detected from the header.
        Returns:
            GaussianOutputFile object
            (optional) the lines of the file
        """
        if not all(isinstance(job, JobType) for job in job_types):
            raise TypeError(f"invalid job type {job}")

        lines = super().read_file(filename)
        header = parse.search_for_block(lines, "#p", "----")

        #### automatically assign job types based on header
        if len(job_types) == 0:
            for name, member in JobType.__members__.items():
                if re.search(f" {member.value}", header):
                    job_types.append(member)

        #### extract parameters
        success = 0
        for line in lines:
            if line.strip().startswith("Normal termination"):
                success += 1

        (geometries, atom_list, energies, scf_iterations,) = parse.read_geometries_and_energies(lines)
        atomic_numbers = list(map(get_number, atom_list))
        bonds = parse.read_bonds(lines)
        charge = int(parse.find_parameter(lines, "Multiplicity", expected_length=6, which_field=2)[0])
        multip = int(parse.find_parameter(lines, "Multiplicity", expected_length=6, which_field=5)[0])

        f = GaussianFile(atomic_numbers, geometries, bonds, job_types=job_types, charge=charge, multiplicity=multip)
        f.energies = energies
        f.scf_iterations = scf_iterations
        f.header = header
        f.success = success

        #### now for some job-type specific attributes
        if JobType.OPT in job_types:
            f.rms_forces = parse.find_parameter(lines, "RMS\s+Force", expected_length=5, which_field=2)
            f.rms_displacements = parse.find_parameter(lines, "RMS\s+Displacement", expected_length=5, which_field=2)

        if JobType.FREQ in job_types:
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

        if return_lines:
            return f, lines
        else:
            return f

    def get_molecule(self, num=None):
        """
        Returns the last molecule (from an optimization job) or the only molecule (from other jobs).

        If ``num`` is specified, returns that job (1-indexed for positive numbers). So ``job.get_molecule(3)`` will return the 3rd element of ``job.molecules``, not the 4th.
        """
        # some methods pass num=None, which overrides setting the default above
        if num is None:
            num = -1

        if not isinstance(num, int):
            raise TypeError("num must be int")

        #### enforce 1-indexing for positive numbers
        if num > 0:
            num += -1

        return self.molecules.molecules[num]
