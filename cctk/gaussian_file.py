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

    All jobs have type ``SP`` by default.
    """

    SP = "sp"
    OPT = "opt"
    FREQ = "freq"
    IRC = "irc"
    NMR = "nmr"
    POP = "pop"
    FORCE = "force"

#### This static variable tells what properties are expected from each JobType.
EXPECTED_PROPERTIES = {
    "sp": ["energy", "scf_iterations",],
    "opt": ["rms_displacement", "rms_force", ],
    "freq": ["gibbs_free_energy", "enthalpy", "frequencies",],
    "nmr": ["isotropic_shielding",],
    "pop": [],
    "force": ["forces",],
}


class GaussianFile(File):
    """
    Class for Gaussian files. Composes ``ConformationalEnsemble``.

    Attributes:
        ensemble (Ensemble): ``ConformationalEnsemble`` instance
        job_types (list): list of `job_type` instances
        route_card (str): optional, route card of .gjf file
        link0 (dict): optional, dictionary of Link 0 commands (e.g. {"mem": "32GB", "nprocshared": 16})
        footer (str): optional, footer of .gjf file
        success (int): number of successful terminations (should be 1 for an opt, 2 for opt and then freq, 1 for a single point energy, etc)
        title (str): optional, title of .gjf file
    """

    def __init__(
        self, job_types, route_card=None, link0=None, footer=None, title="title", success=0
    ):
        """
        Create new GaussianFile object.

	Args:
            job_types (list): list of ``job_type`` instances
            route_card (str): optional, route card of ``.gjf`` file
            link0 (dict): optional, Link 0 commands of ``.gjf`` file
            footer (str): optional, footer of ``.gjf`` file
            title (str): optional, title of ``.gjf`` file
            success (int): num successful terminations
	"""

        if route_card and not isinstance(route_card, str):
            raise TypeError("route card needs to be a string")

        if link0 and not isinstance(link0, dict):
            raise TypeError("link0 needs to be a dict")

        if footer and not isinstance(footer, str):
            raise TypeError("footer needs to be a string")

        if title and not isinstance(title, str):
            raise TypeError("title needs to be a string")

        if success and not isinstance(success, int):
            raise TypeError("success needs to be an integer")

        if job_types is not None:
            if not all(isinstance(job, JobType) for job in job_types):
                raise TypeError(f"invalid job type {job}")

        self.ensemble = ConformationalEnsemble()
        self.route_card = route_card
        self.link0 = link0
        self.footer = footer
        self.title = title
        self.job_types = job_types
        self.success = success

    def __str__(self):
        return f"GaussianFile (title=\"{str(self.title)}\", {len(self.ensemble)} entries in Ensemble)"

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

        assert re.match(r"^#p", route_card), f"route card doesn't start with #p: {route_card}"

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
                text += f"{Z:>2}       {line[0]:>13.8f} {line[1]:>13.8f} {line[2]:>13.8f}\n"
            else:
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
                Default is -1 (e.g. the last molecule), but positive integers will select from self.ensemble(1-indexed).
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
        if (JobType.FREQ in self.job_types) and (self.ensemble[-1:,"frequencies"] is not None):
            freqs = self.ensemble[-1:,"frequencies"]
            if not isinstance(freqs, list) or len(freqs) == 0:
                return list()
            else:
                return list(map(int, np.array(freqs)[np.array(freqs) < 0]))
        else:
            return list()

    @classmethod
    def read_file(cls, filename, return_lines=False):
        """
        Reads a Gaussian``.out`` or ``.gjf`` file and populates the attributes accordingly.
        Only footers from ``opt=modredundant`` can be read automatically --  ``genecep`` custom basis sets, &c must be specified manually.

        Note:

        Will throw ``ValueError`` if there have been no successful iterations.

        Args:
            filename (str): path to the out file
            return_lines (Bool): whether the lines of the file should be returned
        Returns:
            ``GaussianFile`` object (or list of ``GaussianFile`` objects for Link1 files)
            (optional) the lines of the file (or list of lines of file for Link1 files)
        """
        if re.search("gjf$", filename) or re.search("com$", filename):
            return cls._read_gjf_file(filename, return_lines)

        link1_lines = parse.split_link1(super().read_file(filename))
        files = []
        for link1idx, lines in enumerate(link1_lines):
            #### automatically assign job types based on header
            header = parse.search_for_block(lines, "#p", "----")
            if header is None:
                raise ValueError("can't find route card! (perhaps '#p' wasn't employed?)")
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
            charge = parse.find_parameter(lines, "Multiplicity", expected_length=4, which_field=1, split_on="=")[0]
            multip = parse.find_parameter(lines, "Multiplicity", expected_length=4, which_field=3, split_on="=")[0]

            f = GaussianFile(job_types=job_types, route_card=header, link0=link0, footer=footer, success=success)

            molecules = [None] * len(geometries)
            properties = [{} for _ in range(len(geometries))]
            for idx, geom in enumerate(geometries):
                molecules[idx] = Molecule(atomic_numbers, geom, charge=charge, multiplicity=multip, bonds=bonds)
                properties[idx]["energy"] = energies[idx]
                properties[idx]["scf_iterations"] = scf_iterations[idx]
                properties[idx]["link1_idx"] = link1idx
                properties[idx]["filename"] = filename

            #### now for some job-type specific attributes
            if JobType.OPT in job_types:
                rms_forces = parse.find_parameter(lines, "RMS\s+Force", expected_length=5, which_field=2)
                rms_displacements = parse.find_parameter(lines, "RMS\s+Displacement", expected_length=5, which_field=2)

                for idx, force in enumerate(rms_forces):
                    properties[idx]["rms_force"] = force
                    properties[idx]["rms_displacement"] = rms_displacements[idx]

            if JobType.FREQ in job_types:
                enthalpies = parse.find_parameter(lines, "thermal Enthalpies", expected_length=7, which_field=6)
                if len(enthalpies) == 1:
                    properties[-1]["enthalpy"] = enthalpies[0]
                elif len(enthalpies) > 1:
                    raise ValueError("too many enthalpies found!")

                gibbs_vals = parse.find_parameter(lines, "thermal Free Energies", expected_length=8, which_field=7)
                if len(gibbs_vals) == 1:
                    properties[-1]["gibbs_free_energy"] = gibbs_vals[0]
                elif len(gibbs_vals) > 1:
                    raise ValueError("too many gibbs free energies found!")

                frequencies = []
                try:
                    frequencies += parse.find_parameter(lines, "Frequencies", expected_length=5, which_field=2)
                    frequencies += parse.find_parameter(lines, "Frequencies", expected_length=5, which_field=3)
                    frequencies += parse.find_parameter(lines, "Frequencies", expected_length=5, which_field=4)
                    properties[-1]["frequencies"] = sorted(frequencies)
                except:
                    raise ValueError("error finding frequencies")

            if JobType.NMR in job_types:
                assert len(molecules) == 1, "NMR jobs should not be combined with optimizations!"
                nmr_shifts = parse.read_nmr_shifts(lines, molecules[0].num_atoms())
                properties[0]["isotropic_shielding"] = nmr_shifts.view(OneIndexedArray)

            if JobType.FORCE in job_types:
                assert len(molecules) == 1, "force jobs should not be combined with optimizations!"
                forces = parse.read_forces(lines)
                properties[0]["forces"] = forces

            if JobType.POP in job_types:
                if re.search("hirshfeld", f.route_card):
                    charges, spins = parse.read_hirshfeld_charges(lines)
                    properties[-1]["hirshfeld_charges"] = charges
                    properties[-1]["hirshfeld_spins"] = spins

            try:
                charges = parse.read_mulliken_charges(lines)
                properties[-1]["mulliken_charges"] = charges
            except:
                pass

            try:
                dipole = parse.read_dipole_moment(lines)
                properties[-1]["dipole_moment"] = dipole
            except:
                pass

            for mol, prop in zip(molecules, properties):
                f.ensemble.add_molecule(mol, properties=prop)

            f.check_has_properties()
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
        Reads a Gaussian ``.gjf`` or ``.com`` file and populates the attributes accordingly.

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
        geometry = []

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
                    geometry.append([pieces[1], pieces[2], pieces[3]])

            if (in_geom == False) and (len(geometry) > 0):
                if footer:
                    footer = footer + "\n" + line
                else:
                    if len(line.strip()) > 0:
                        footer = line

        try:
            atomic_numbers = np.array(atomic_numbers, dtype=np.int8)
        except:
            atomic_numbers = np.array(list(map(get_number, atomic_numbers)), dtype=np.int8)

        job_types = cls._assign_job_types(header)

        f = GaussianFile(job_types=job_types, route_card=header, link0=link0, footer=footer, title=title)
        f.ensemble.add_molecule(Molecule(atomic_numbers, geometry, charge=charge, multiplicity=multip))
        if return_lines:
            return f, lines
        else:
            return f

    def get_molecule(self, num=None):
        """
        Returns the last molecule (from an optimization job) or the only molecule (from other jobs).

        If ``num`` is specified, returns ``self.ensemble.molecule_list()[num]``
        """
        # some methods pass num=None, which overrides setting the default above
        if num is None:
            num = -1

        if not isinstance(num, int):
            raise TypeError("num must be int")

        return self.ensemble.molecule_list()[num]

    @classmethod
    def _assign_job_types(cls, header):
        """
        Assigns ``JobType`` objects from route card. ``Job.Type.SP`` is assigned by default.

        For instance, "#p opt freq=noraman" would give an output of ``[JobType.SP, JobType.OPT, JobType.FREQ]``.

        Args:
            header (str): Gaussian route card

        Returns:
            list of ``JobType`` objects
        """
        job_types = []
        for name, member in JobType.__members__.items():
            if re.search(f" {member.value}", str(header), re.IGNORECASE):
                job_types.append(member)
        if JobType.SP not in job_types:
            job_types.append(JobType.SP)
        return job_types

    def check_has_properties(self):
        """
        Checks that the file has all the appropriate properties for its job types, and raises ValueError if not.

        This only checks the last molecule in ``self.ensemble``, for now.
        """
        return True
#        for job_type in self.job_types:
#            for prop in EXPECTED_PROPERTIES[job_type.value]:
#                if not self.ensemble.has_property(-1, prop):
#                    raise ValueError(f"expected property {prop} for job type {job_type}, but it's not there!")

    @classmethod
    def write_ensemble_to_file(cls, filename, ensemble, route_card, link0={"mem": "32GB", "nprocshared": 16}, footer=None, title="title", print_symbol=False):
            """
            Write each structure in the specified ensemble to a single Gaussian input file
            by using the Link1 specification.

            Args:
                filename (str): where to write the file
                ensemble (Ensemble): ``Ensemble`` object to write
                route_card (str or list): to use the same route card for every link, use a single string;
                                          otherwise, provide a list whose entries parallel the ensemble members
                link0 (dict or list of dicts): to use the same memory/processors for every link, use a single string;
                                               otherwise, provide a list
                footer (None/str or list): use None for no text after geometry, provide a str to specify a footer,
                                           or provide some combination of the above as a list
                title (str or list): use a single string to provide a generic title for every link or a list as above
                print_symbol (bool or list): whether to print atomic symbols or atomic numbers in the geometry specification;
                                             use a single bool or a list as above

            """
            n_geometries = len(ensemble)
            assert len(ensemble) > 0, "cannot write a blank ensemble"

            if isinstance(route_card, str):
                assert re.match(r"^#p", route_card), "route card should start with #p: {route_card}"
                route_card = [route_card for m in ensemble._items]
            elif isinstance(route_card, list):
                assert len(route_card) == n_geometries, f"expected {n_geometries} route cards but got {len(route_card)}"
                for card in route_card:
                    assert isinstance(card, str), "expected route card to be a str"
                    assert re.match(r"^#p", route_card), "route card should start with #p: f{card}"
            else:
                raise ValueError(f"unexpected type for route_card: {str(type(route_card))}")

            if isinstance(link0, dict):
                link0 = [link0 for m in ensemble._items]
            elif isinstance(link0, list):
                assert len(link0) == n_geometries, f"expected {n_geometries} link0 entries, but got {len(link0)}"
                for d in link0:
                    assert isinstance(d, dict), f"expected dict for link0 but got {str(type(d))}"
            else:
                raise ValueError(f"unexpected type for link0: {str(type(link0))}")

            if footer is None or isinstance(footer, str):
                footer = [footer for m in ensemble._items]
            elif isinstance(footer, list):
                assert len(footer) == n_geometries, f"expected {n_geometries} footers, but got {len(footer)}"
                for f in footer:
                    assert f is None or isinstance(f, str), f"expected str or None for footer but got {str(type(f))}"
            else:
                raise ValueError(f"unexpected type for footer: {str(type(footer))}")

            if isinstance(title, str):
                assert len(title.strip()) > 0, "zero-length titles not allowed"
                title = [title for m in ensemble._items]
            elif isinstance(title, list):
                assert len(title) == n_geometries, f"expected {n_geometries} route cards but got {len(title)}"
                for card in title:
                    assert isinstance(card, str), "expected title to be a str"
                    assert len(title.strip()) > 0, "zero-length titles are not allowed"
            else:
                raise ValueError(f"unexpected type for title: {str(type(title))}")

            if isinstance(print_symbol, bool):
                print_symbol = [print_symbol for m in ensemble._items]
            elif isinstance(print_symbol, list):
                assert len(print_symbol) == n_geometries, f"expected {n_geometries} print_symbol entries but got {len(print_symbol)}"
                for s in print_symbol:
                    assert isinstance(s, bool), f"expected bool for print_symbol but got {str(type(s))}"
            else:
                raise ValueError(f"unexpected type for print_symbol: {str(type(print_symbol))}")

            for idx, molecule in enumerate(ensemble._items):
                if idx == 0:
                    cls.write_molecule_to_file(filename, molecule, route_card[idx], link0[idx], footer=footer[idx], title=title[idx], print_symbol=print_symbol[idx], append=False)
                else:
                    cls.write_molecule_to_file(filename, molecule, route_card[idx], link0[idx], footer=footer[idx], title=title[idx], print_symbol=print_symbol[idx], append=True)


