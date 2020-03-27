import sys, re
import numpy as np

from enum import Enum

from cctk import File, Ensemble, Molecule, ConformationalEnsemble, OneIndexedArray
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
    IRC = "irc"
    NMR = "nmr"
    POP = "pop"
    FORCE = "force"


class GaussianFile(File):
    """
    Class for Gaussian files. Composes ``ConformationalEnsemble``.

    Attributes:
        molecules (Ensemble): ``ConformationalEnsemble`` instance
        job_types (list): list of `job_type` instances
        route_card (str): optional, route card of .gjf file
        link0 (dict): optional, dictionary of Link 0 commands (e.g. {"mem": "32GB", "nprocshared": 16})
        footer (str): optional, footer of .gjf file
        success (int): number of successful terminations (should be 1 for an opt, 2 for opt and then freq, 1 for a single point energy, etc)
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

    def __init__(
        self, atomic_numbers, geometries, bonds=None, charge=None, multiplicity=None, job_types=None, route_card=None, link0=None, footer=None, title="title",
    ):
        """
        Create new GaussianFile object.

		Args:
            atomic_numbers (list): list of atomic numbers
            geometries (list): list of lists of 3-tuples of xyz coordinates
            bonds (nx.Graph): Graph object containing connectivity information (1-indexed)
            charges (int): the charge of the molecules
            multiplicities (int): the spin states of the molecules (1 corresponds to singlet, 2 to doublet, 3 to triplet, etc. -- so a multiplicity of 1 is equivalent to S=0)
            job_types (list): list of ``job_type`` instances
            route_card (str): optional, route card of ``.gjf`` file
            link0 (dict): optional, Link 0 commands of ``.gjf`` file
            footer (str): optional, footer of ``.gjf`` file
            title (str): optional, title of ``.gjf`` file
		"""

        if route_card and not isinstance(route_card, str):
            raise TypeError("route card needs to be a string")

        if link0 and not isinstance(link0, dict):
            raise TypeError("link0 needs to be a dict")

        if footer and not isinstance(footer, str):
            raise TypeError("footer needs to be a string")

        if title and not isinstance(title, str):
            raise TypeError("title needs to be a string")

        if job_types is not None:
            if not all(isinstance(job, JobType) for job in job_types):
                raise TypeError(f"invalid job type {job}")

        if (len(atomic_numbers) > 0) and (len(geometries) > 0):
            arguments = {"atomic_numbers": atomic_numbers, "geometries": geometries, "bonds": bonds, "charge": charge, "multiplicity": multiplicity}
            self.molecules = ConformationalEnsemble(**arguments)

        self.route_card = route_card
        self.link0 = link0
        self.footer = footer
        self.title = title
        self.job_types = job_types

    @classmethod
    def write_molecule_to_file(cls, filename, molecule, route_card, link0={"mem": "32GB", "nprocshared": 16}, footer=None, title="title", append=False, print_symbol=False):
        """
        Write a ``.gjf`` file using the given molecule.

        Args:
            filename (str): path to the new file
            molecule (Molecule): which molecule to use -- a ``Molecule`` object.
            route_card (str): route card for new file
            link0 (dict): dictionary of Link 0 commands
            footer (str): footer for new file
            title (str): title of the file, defaults to "title"
            append (Bool): whether or not to append to file using Link1 specifications
            print_symbol (Bool): whether to print atomic symbols (instead of atomic numbers)
        """
        if not isinstance(molecule, Molecule):
            raise TypeError("need a valid molecule to write a file!")

        if (route_card is None) or (not isinstance(route_card, str)):
            raise ValueError("can't write a file without a route card")

        #### generate the text
        text = ""
        if append:
            text += "--Link1--\n"

        if isinstance(link0, dict):
            for key, val in link0.items():
                text += f"%{key}={val}\n"

        text += f"{route_card.strip()}\n\n{title}\n\n"

        text += f"{int(molecule.charge)} {int(molecule.multiplicity)}\n"
        for index, Z in enumerate(molecule.atomic_numbers, start=1):
            line = molecule.get_vector(index)
            if print_symbol:
                Z = get_symbol(Z)
            text += f"{Z:2d}       {line[0]:>13.8f} {line[1]:>13.8f} {line[2]:>13.8f}\n"

        text += "\n"
        if footer is not None:
            text += f"{footer.strip()}\n\n"

        #### write the file
        if append:
            super().append_to_file(filename, text)
        else:
            super().write_file(filename, text)

    def write_file(self, filename, molecule=None, route_card=None, link0=None, footer=None, **kwargs):
        """
        Write a ``.gjf`` file, using object attributes. If no header/footer is specified, the object's header/footer will be used.

        Args:
            filename (str): path to the new file
            molecule (int): which molecule to use -- passed to ``self.get_molecule()``.
                Default is -1 (e.g. the last molecule), but positive integers will select from self.molecules (1-indexed).
                A ``Molecule`` object can also be passed, in which case that molecule will be written to the file.
            route_card (str): route card for new file
            link0 (dict): dictionary of Link 0 commands (e.g. {"mem": "32GB", "nprocshared": 16}
            footer (str): footer for new file
        """
        if not isinstance(molecule, Molecule):
            molecule = self.get_molecule(molecule)

        if route_card is None:
            route_card = self.route_card

        if link0 is None:
            link0 = self.link0

        if footer is None:
            footer = self.footer

        self.write_molecule_to_file(filename, molecule, route_card, link0, footer, **kwargs)

    def num_imaginaries(self):
        """
        Returns the number of imaginary frequencies.
        """
        return len(self.imaginaries())

    def imaginaries(self):
        """
        Returns the imaginary frequencies, rounded to the nearest integer.
        """
        if JobType.FREQ in self.job_types:
            return list(map(int, np.array(self.frequencies)[np.array(self.frequencies) < 0]))
        else:
            raise TypeError("not a frequency job! can't get # imaginary frequencies!")

    @classmethod
    def read_file(cls, filename, return_lines=False):
        """
        Reads a Gaussian``.out`` or ``.gjf`` file and populates the attributes accordingly.
        Only footers from ``opt=modredundant`` can be read automatically --  ``genecep`` custom basis sets, &c must be specified manually.

        Will throw ``ValueError`` if there have been no successful iterations.

        Args:
            filename (str): path to the out file
            return_lines (Bool): whether the lines of the file should be returned
        Returns:
            ``GaussianFile`` object (or list of ``GaussianFile`` objects for Link1 files)
            (optional) the lines of the file (or list of lines of file for Link1 files)
        """
        if re.search("gjf$", filename):
            return cls._read_gjf_file(filename, return_lines)

        link1_lines = parse.split_link1(super().read_file(filename))
        files = []
        for lines in link1_lines:
            #### automatically assign job types based on header
            header = parse.search_for_block(lines, "#p", "----")
            job_types = cls._assign_job_types(header)

            link0 = parse.extract_link0(lines)

            #### extract parameters
            success = 0
            for line in lines:
                if line.strip().startswith("Normal termination"):
                    success += 1

            (geometries, atom_list, energies, scf_iterations,) = parse.read_geometries_and_energies(lines)
            atomic_numbers = []

            #### convert to right datatype
            try:
                atomic_numbers = np.array(atom_list, dtype=np.int8)
            except:
                atomic_numbers = np.array(list(map(get_number, atom_list)), dtype=np.int8)

            footer = None
            if re.search("modredundant", str(header)):
                footer = parse.search_for_block(lines, "^ The following ModRedundant input section", "^ $", count=1, join="\n")
                footer = "\n".join(list(footer.split("\n"))[1:])  # get rid of the first line
                footer = "\n".join([" ".join(list(filter(None, line.split(" ")))) for line in footer.split("\n")])

            bonds = parse.read_bonds(lines)
            charge = parse.find_parameter(lines, "Multiplicity", expected_length=6, which_field=2)[0]
            multip = parse.find_parameter(lines, "Multiplicity", expected_length=6, which_field=5)[0]

            f = GaussianFile(atomic_numbers, geometries, bonds, job_types=job_types, charge=charge, multiplicity=multip)
            f.energies = energies
            f.scf_iterations = scf_iterations
            f.route_card = header
            f.link0 = link0
            f.footer = footer
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

            if JobType.NMR in job_types:
                assert len(f.molecules) == 1, "NMR jobs should not be combined with optimizations!"
                nmr_shifts = parse.read_nmr_shifts(lines, f.molecules[-1].num_atoms())
                f.molecules[0].nmr_isotropic = nmr_shifts.view(OneIndexedArray)

            files.append(f)

        if return_lines:
            if len(link1_lines) == 1:
                return files[0], link1_lines[0]
            else:
                return files, link1_lines
        else:
            if len(link1_lines) == 1:
                return files[0]
            else:
                return files

    @classmethod
    def _read_gjf_file(cls, filename, return_lines=False):
        """
        Reads a Gaussian ``.gjf`` file and populates the attributes accordingly.

        Args:
            filename (str): path to the out file
            return_lines (Bool): whether the lines of the file should be returned
        Returns:
            GaussianFile object
            (optional) the lines of the file
        """
        lines = super().read_file(filename)
        header = None
        link0 = {}
        footer = None
        header_done = False
        title = None
        charge = None
        multip = None
        in_geom = False
        atomic_numbers = []
        geometries = []

        for idx, line in enumerate(lines):
            if header is None:
                if re.match("\%", line):
                    pieces = line[1:].split("=")
                    link0[pieces[0]] = pieces[1]
                    continue
                if re.match("#", line):
                    header = line
                    continue

            if (title is None) and (header is not None):
                if header_done:
                    if len(line.strip()) > 0:
                        title = line
                else:
                    if len(line.strip()) > 0:
                        header = header + line
                    else:
                        header_done = True
                continue

            if (title is not None) and (charge is None):
                if len(line.strip()) > 0:
                    pieces = list(filter(None, line.split(" ")))
                    assert len(pieces) == 2, f"can't parse line {line}"

                    charge = int(pieces[0])
                    multip = int(pieces[1])
                    in_geom = True
                    continue

            if in_geom == True:
                if len(line.strip()) == 0:
                    in_geom = False
                else:
                    pieces = list(filter(None, line.split(" ")))
                    assert len(pieces) == 4, f"can't parse line {line}"

                    atomic_numbers.append(pieces[0])
                    geometries.append([pieces[1], pieces[2], pieces[3]])

            if (in_geom == False) and (len(geometries) > 0):
                if footer:
                    footer = footer + "\n" + line
                else:
                    if len(line.strip()) > 0:
                        footer = line

        try:
            atomic_numbers = np.array(atomic_numbers, dtype=np.int8)
        except:
            atomic_numbers = np.array(list(map(get_number, atomic_numbers)), dtype=np.int8)

        geometries = np.array([geometries])
        job_types = cls._assign_job_types(header)

        f = GaussianFile(atomic_numbers, geometries, job_types=job_types, charge=charge, multiplicity=multip)
        f.route_card = header
        f.link0 = link0
        f.footer = footer
        f.title = title
        f.success = 0

        if return_lines:
            return f, lines
        else:
            return f

    def get_molecule(self, num=None):
        """
        Returns the last molecule (from an optimization job) or the only molecule (from other jobs).

        If ``num`` is specified, returns ``self.molecules[num]``
        """
        # some methods pass num=None, which overrides setting the default above
        if num is None:
            num = -1

        if not isinstance(num, int):
            raise TypeError("num must be int")

        return self.molecules[num]

    @classmethod
    def write_ensemble_to_file(cls, filename, ensemble, route_cards, kwargs):
        """
        Writes an Ensemble to a file using Link1 specification.

        Args:
            filename (str): where to write the file
            ensemble (Ensemble): ``Ensemble`` object to write
            headers (list): headers for each ``write_molecule_to_file`` call
            kwargs (list of dict): arguments for each ``write_molecule_to_file`` call
        """
        for idx, molecule in enumerate(ensemble.molecules):
            if idx == 0:
                cls.write_molecule_to_file(filename, molecule, route_cards[idx], append=False, **kwargs[idx])
            else:
                cls.write_molecule_to_file(filename, molecule, route_cards[idx], append=True, **kwargs[idx])


    @classmethod
    def _assign_job_types(cls, header):
        """
        Assigns ``JobType`` objects from route card.

        For instance, "#p opt freq=noraman" would give an output of [JobType.OPT, JobType.FREQ].

        Args:
            header (str): Gaussian route card

        Returns:
            list of ``JobType`` objects
        """
        job_types = []
        for name, member in JobType.__members__.items():
            if re.search(f" {member.value}", str(header), re.IGNORECASE):
                job_types.append(member)
        return job_types
