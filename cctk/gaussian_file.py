import re, warnings
import numpy as np

from enum import Enum

from cctk import File, Molecule, ConformationalEnsemble, OneIndexedArray
from cctk.helper_functions import get_symbol, get_number, get_corrected_free_energy
import cctk

import cctk.parse_gaussian as parse


class GaussianJobType(Enum):
    """
    Class representing allowed Gaussian job types. Not an exhaustive list, but should be fairly comprehensive.

    The value should be the Gaussian keyword, to permit automatic assignment.

    All jobs have type ``SP`` by default.
    """

    SP = "sp"
    """
    Single point energy calculation.
    """

    OPT = "opt"
    """
    Geometry optimization.
    """

    FREQ = "freq"
    """
    Hessian calculation.
    """

    IRC = "irc"
    """
    Intrinsic reaction coordinate calculation.
    """

    NMR = "nmr"
    """
    NMR shielding prediction.
    """

    POP = "pop"
    """
    Population analysis.
    """

    FORCE = "force"
    """
    Gradient calculation.
    """

#### This static variable tells what properties are expected from each JobType.
EXPECTED_PROPERTIES = {
    "sp": ["energy", "scf_iterations",],
    "opt": ["rms_displacement", "rms_force",],
    "freq": ["gibbs_free_energy", "enthalpy", "frequencies",],
    "nmr": ["isotropic_shielding",],
    "pop": [],
    "force": ["forces",],
}


class GaussianFile(File):
    """
    Class representing Gaussian input/output files.

    Attributes:
        ensemble (ConformationalEnsemble): ``ConformationalEnsemble`` instance
        job_types (list): list of `job_type`` instances
        route_card (str): optional, route card of .gjf file
        link0 (dict): optional, dictionary of Link 0 commands (e.g. {"mem": "32GB", "nprocshared": 16})
        footer (str): optional, footer of .gjf file
        successful_terminations (int): number of successful terminations (should be 1 for an opt, 2 for opt and then freq, 1 for a single point energy, etc)
        elapsed_time (float): total time for job in seconds
        title (str): optional, title of .gjf file
    """

    def __init__(
        self, job_types=None, route_card=None, link0=None, footer=None, title="title", success=0, elapsed_time=0.0, molecule=None,
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
            elapsed_time (float): total time for job in seconds
            molecule (cctk.Molecule): molecule to initiate, if desired
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

        if not isinstance(elapsed_time, (float, int)) or elapsed_time < 0.0:
            raise TypeError(f"elapsed_time invalid: {elapsed_time}")

        if job_types is not None:
            if isinstance(job_types, str):
                raise ValueError(f"invalid job_types {job_types} - did you mean to call GaussianFile.read_file({job_types})?")
            if not all(isinstance(job, GaussianJobType) for job in job_types):
                raise TypeError(f"invalid job_types {job_types}")

        self.ensemble = ConformationalEnsemble()

        if molecule is not None:
            assert isinstance(molecule, Molecule), "molecule is not a valid cctk.Molecule!"
            self.ensemble.add_molecule(molecule)

        self.route_card = route_card
        self.link0 = link0
        self.footer = footer
        self.title = title
        self.job_types = job_types
        self.successful_terminations = success
        self.elapsed_time = elapsed_time

    def __str__(self):
        return f"GaussianFile (title=\"{str(self.title)}\", {len(self.ensemble)} entries in Ensemble)"

    @classmethod
    def write_molecule_to_file(cls, filename, molecule, route_card, link0={"mem": "32GB", "nprocshared": 16}, footer=None, title="title", append=False, print_symbol=False, point_charges=None):
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

        if not re.match(r"^#p", route_card):
            warnings.warn(f"route card doesn't start with #p: {route_card}")

        if point_charges is not None:
            assert isinstance(point_charges, list), "point_charges must be list"
            assert all([isinstance(pc, cctk.PointCharge) for pc in point_charges]), "point_charges must be list of point charges"
            assert re.search(r"charge", route_card, flags=re.IGNORECASE), "charge must be in route_card if point_charges are present"

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

        if point_charges is not None:
            for point_charge in point_charges:
                text += f"{point_charge.coordinates[0]:>13.8f} {point_charge.coordinates[1]:>13.8f} {point_charge.coordinates[2]:>13.8f} {point_charge.charge:.5f}\n"
            text += "\n"

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
        if (GaussianJobType.FREQ in self.job_types) and (self.ensemble[-1:,"frequencies"] is not None):
            freqs = self.ensemble[-1:,"frequencies"]
            if not isinstance(freqs, list) or len(freqs) == 0:
                return list()
            else:
                return list(map(int, np.array(freqs)[np.array(freqs) < 0]))
        else:
            return list()

    @classmethod
#    @profile
    def read_file(cls, filename, return_lines=False, extended_opt_info=False):
        """
        Reads a Gaussian``.out`` or ``.gjf`` file and populates the attributes accordingly.
        Only footers from ``opt=modredundant`` can be read automatically --  ``genecep`` custom basis sets, &c must be specified manually.

        Note:

        Will throw ``ValueError`` if there have been no successful iterations.

        Args:
            filename (str): path to the out file
            return_lines (Bool): whether the lines of the file should be returned
            extended_opt_info (Bool): if full parameters about each opt step should be collected
                (by default, only ``rms_displacement`` and ``rms_force`` are collected)
        Returns:
            ``GaussianFile`` object (or list of ``GaussianFile`` objects for Link1 files)
            (optional) the lines of the file (or list of lines of file for Link1 files)
        """
        if re.search("gjf$", filename) or re.search("com$", filename):
            return cls._read_gjf_file(filename, return_lines)

        link1_lines = parse.split_link1(filename)
        files = []

        for link1idx, lines in enumerate(link1_lines):
            #### automatically assign job types based on header
            header = lines.search_for_block("#p", "----", format_line=lambda x: x.lstrip(), join="")
            if header is None:
                raise ValueError("can't find route card! (perhaps '#p' wasn't employed?)")
            job_types = cls._assign_job_types(header)

            link0 = parse.extract_link0(lines)

            title = ""
            title_block = lines.search_for_block("l101.exe", "Symbolic Z-matrix", join="\n")
            if title_block is not None:
                for line in title_block.split("\n")[1:]:
                    if not re.search("-----", line):
                        title += line


            (geometries, atom_list, energies, scf_iterations, success, elapsed_time) = parse.read_geometries_and_energies(lines)
            success, elapsed_time = parse.extract_success_and_time(lines)
            atomic_numbers = []

            #### convert to right datatype
            try:
                atomic_numbers = np.array(atom_list, dtype=np.int8)
            except Exception as e:
                atomic_numbers = np.array(list(map(get_number, atom_list)), dtype=np.int8)

            footer = None
            if re.search("modredundant", str(header)):
                footer = lines.search_for_block("^ The following ModRedundant input section", "^ $", count=1, join="\n")
                if footer is not None:
                    footer = "\n".join(list(footer.split("\n"))[1:])  # get rid of the first line
                    footer = "\n".join([" ".join(list(filter(None, line.split(" ")))) for line in footer.split("\n")])

            bonds = parse.read_bonds(lines)
            charge, multip =  lines.find_parameter("Multiplicity", expected_length=4, which_field=[1,3], split_on="=")[0]

            f = GaussianFile(job_types=job_types, route_card=header, link0=link0, footer=footer, success=success, elapsed_time=elapsed_time, title=title)

            molecules = [None] * len(geometries)
            properties = [{} for _ in range(len(geometries))]
            for idx, geom in enumerate(geometries):
                molecules[idx] = Molecule(atomic_numbers, geom, charge=charge, multiplicity=multip, bonds=bonds)
                if idx < len(energies):
                    properties[idx]["energy"] = energies[idx]
                if idx < len(scf_iterations):
                    properties[idx]["scf_iterations"] = scf_iterations[idx]
                properties[idx]["link1_idx"] = link1idx
                properties[idx]["filename"] = filename
                properties[idx]["iteration"] = idx

            #### now for some job-type specific attributes
            if GaussianJobType.OPT in job_types:
                rms_forces = lines.find_parameter("RMS\s+Force", expected_length=5, which_field=2)
                rms_displacements = lines.find_parameter("RMS\s+Displacement", expected_length=5, which_field=2)

                if extended_opt_info:
                    max_forces = lines.find_parameter("Maximum Force", expected_length=5, which_field=2)
                    max_displacements = lines.find_parameter("Maximum Displacement", expected_length=5, which_field=2)
                    max_gradients = lines.find_parameter("Cartesian Forces:", expected_length=6, which_field=3)
                    rms_gradients = lines.find_parameter("Cartesian Forces:", expected_length=6, which_field=5)
                    max_int_forces = lines.find_parameter("Internal  Forces:", expected_length=6, which_field=3)
                    rms_int_forces = lines.find_parameter("Internal  Forces:", expected_length=6, which_field=5)
                    delta_energy = lines.find_parameter("Predicted change in Energy", expected_length=4, which_field=3, cast_to_float=False)

                for idx, force in enumerate(rms_forces):
                    properties[idx]["rms_force"] = force
                    properties[idx]["rms_displacement"] = rms_displacements[idx]

                    if extended_opt_info:
                        if idx < len(max_forces):
                            properties[idx]["max_force"] = max_forces[idx]

                        if idx < len(max_displacements):
                            properties[idx]["max_displacement"] = max_displacements[idx]

                        if idx < len(max_gradients):
                            properties[idx]["max_gradient"] = max_gradients[idx]

                        if idx < len(rms_gradients):
                            properties[idx]["rms_gradient"] = rms_gradients[idx]

                        if idx < len(max_int_forces):
                            properties[idx]["max_internal_force"] = max_int_forces[idx]

                        if idx < len(rms_int_forces):
                            properties[idx]["rms_internal_force"] = rms_int_forces[idx]

                        if idx < len(delta_energy):
                            change_in_energy = re.sub(r"Energy=", "", delta_energy[idx])
                            properties[idx]["predicted_change_in_energy"] = float(change_in_energy.replace('D', 'E'))

            if GaussianJobType.FREQ in job_types:
                enthalpies = lines.find_parameter("thermal Enthalpies", expected_length=7, which_field=6)
                if len(enthalpies) == 1:
                    properties[-1]["enthalpy"] = enthalpies[0]
                elif len(enthalpies) > 1:
                    raise ValueError(f"unexpected # of enthalpies found!\nenthalpies = {enthalpies}")

                gibbs_vals = lines.find_parameter("thermal Free Energies", expected_length=8, which_field=7)
                if len(gibbs_vals) == 1:
                    properties[-1]["gibbs_free_energy"] = gibbs_vals[0]
                elif len(gibbs_vals) > 1:
                    raise ValueError(f"unexpected # gibbs free energies found!\ngibbs free energies = {gibbs_vals}")

            if GaussianJobType.FREQ in job_types:
                enthalpies = lines.find_parameter("thermal Enthalpies", expected_length=7, which_field=6)
                if len(enthalpies) == 1:
                    properties[-1]["enthalpy"] = enthalpies[0]
                elif len(enthalpies) > 1:
                    raise ValueError(f"unexpected # of enthalpies found!\nenthalpies = {enthalpies}")

                gibbs_vals = lines.find_parameter("thermal Free Energies", expected_length=8, which_field=7)
                if len(gibbs_vals) == 1:
                    properties[-1]["gibbs_free_energy"] = gibbs_vals[0]
                elif len(gibbs_vals) > 1:
                    raise ValueError(f"unexpected # gibbs free energies found!\ngibbs free energies = {gibbs_vals}")

                frequencies = []
                try:
                    frequencies = sum(lines.find_parameter("Frequencies", expected_length=5, which_field=[2,3,4]), [])
                    properties[-1]["frequencies"] = sorted(frequencies)
                except Exception as e:
                    raise ValueError("error finding frequencies")

                #  Temperature   298.150 Kelvin.  Pressure   1.00000 Atm.
                temperature = lines.find_parameter("Temperature", expected_length=6, which_field=1)
                if len(temperature) == 1:
                    properties[-1]["temperature"] = temperature[0]
                    try:
                        corrected_free_energy = get_corrected_free_energy(gibbs_vals[0], frequencies, frequency_cutoff=100.0, temperature=temperature[0])
                        properties[-1]["quasiharmonic_gibbs_free_energy"] = float(f"{float(corrected_free_energy):.6f}") # yes this is dumb
                    except Exception as e:
                        pass


            if GaussianJobType.NMR in job_types:
                nmr_shifts = parse.read_nmr_shifts(lines, molecules[0].num_atoms())
                if nmr_shifts is not None:
                    properties[-1]["isotropic_shielding"] = nmr_shifts.view(OneIndexedArray)

                if re.search("nmr=mixed", f.route_card, flags=re.IGNORECASE) or re.search("nmr=spinspin", f.route_card,flags=re.IGNORECASE):
                    couplings = parse.read_j_couplings(lines, molecules[0].num_atoms())
                    if couplings is not None:
                        properties[-1]["j_couplings"] = couplings

            if GaussianJobType.FORCE in job_types:
                assert len(molecules) == 1, "force jobs should not be combined with optimizations!"
                forces = parse.read_forces(lines)
                properties[0]["forces"] = forces

            if GaussianJobType.POP in job_types:
                if re.search("hirshfeld", f.route_card) or re.search("cm5", f.route_card):
                    charges, spins = parse.read_hirshfeld_charges(lines)
                    properties[-1]["hirshfeld_charges"] = charges
                    properties[-1]["hirshfeld_spins"] = spins

            try:
                charges = parse.read_mulliken_charges(lines)
                properties[-1]["mulliken_charges"] = charges
            except Exception as e:
                pass

            try:
                dipole = parse.read_dipole_moment(lines)
                properties[-1]["dipole_moment"] = dipole
            except Exception as e:
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
        except Exception as e:
            atomic_numbers = np.array(list(map(get_number, atomic_numbers)), dtype=np.int8)

        job_types = cls._assign_job_types(header)

        f = GaussianFile(job_types=job_types, route_card=header, link0=link0, footer=footer, title=title)
        f.ensemble.add_molecule(Molecule(atomic_numbers, geometry, charge=charge, multiplicity=multip))
        if return_lines:
            return f, lines
        else:
            return f

    def get_molecule(self, num=None, properties=False):
        """
        Returns the last molecule (from an optimization job) or the only molecule (from other jobs).

        If ``num`` is specified, returns ``self.ensemble.molecule_list()[num]``
        If ``properties`` is True, returns ``(molecule, properties)``.
        """
        # some methods pass num=None, which overrides setting the default above
        if num is None:
            num = -1
        assert isinstance(num, int), "num must be int"

        if properties:
            return self.ensemble.molecule_list()[num], self.ensemble.properties_list()[num]
        else:
            return self.ensemble.molecule_list()[num]

    @classmethod
    def _assign_job_types(cls, header):
        """
        Assigns ``GaussianJobType`` objects from route card. ``GaussianJobType.SP`` is assigned by default.

        For instance, "#p opt freq=noraman" would give an output of ``[GaussianJobType.SP, GaussianJobType.OPT, GaussianJobType.FREQ]``.

        Args:
            header (str): Gaussian route card

        Returns:
            list of ``GaussianJobType`` objects
        """
        job_types = []
        for name, member in GaussianJobType.__members__.items():
            if re.search(f" {member.value}", str(header), re.IGNORECASE):
                job_types.append(member)
        if GaussianJobType.SP not in job_types:
            job_types.append(GaussianJobType.SP)
        return job_types

    def check_has_properties(self):
        """
        Checks that the file has all the appropriate properties for its job types, and raises ValueError if not.

        This only checks the last molecule in ``self.ensemble``, for now.
        """
        if self.successful_terminations > 0:
            if self.successful_terminations == 1 and ((GaussianJobType.OPT in self.job_types) and (GaussianJobType.FREQ in self.job_types)):
                return # opt freq jobs should have two terminations
            for job_type in self.job_types:
                for prop in EXPECTED_PROPERTIES[job_type.value]:
                    if not self.ensemble.has_property(-1, prop):
                        raise ValueError(f"expected property {prop} for job type {job_type}, but it's not there!")
        else:
            return

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
                route_card = [route_card for _ in ensemble._items]
            elif isinstance(route_card, list):
                assert len(route_card) == n_geometries, f"expected {n_geometries} route cards but got {len(route_card)}"
                for card in route_card:
                    assert isinstance(card, str), "expected route card to be a str"
            else:
                raise ValueError(f"unexpected type for route_card: {str(type(route_card))}")

            if isinstance(link0, dict):
                link0 = [link0 for _ in ensemble._items]
            elif isinstance(link0, list):
                assert len(link0) == n_geometries, f"expected {n_geometries} link0 entries, but got {len(link0)}"
                for d in link0:
                    assert isinstance(d, dict), f"expected dict for link0 but got {str(type(d))}"
            else:
                raise ValueError(f"unexpected type for link0: {str(type(link0))}")

            if footer is None or isinstance(footer, str):
                footer = [footer for _ in ensemble._items]
            elif isinstance(footer, list):
                assert len(footer) == n_geometries, f"expected {n_geometries} footers, but got {len(footer)}"
                for f in footer:
                    assert f is None or isinstance(f, str), f"expected str or None for footer but got {str(type(f))}"
            else:
                raise ValueError(f"unexpected type for footer: {str(type(footer))}")

            if isinstance(title, str):
                assert len(title.strip()) > 0, "zero-length titles not allowed"
                title = [title for _ in ensemble._items]
            elif isinstance(title, list):
                assert len(title) == n_geometries, f"expected {n_geometries} route cards but got {len(title)}"
                for card in title:
                    assert isinstance(card, str), "expected title to be a str"
                    assert len(title.strip()) > 0, "zero-length titles are not allowed"
            else:
                raise ValueError(f"unexpected type for title: {str(type(title))}")

            if isinstance(print_symbol, bool):
                print_symbol = [print_symbol for _ in ensemble._items]
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

    def add_custom_basis_set(self, name, add_all_elements=False, return_string=False):
        """
        Appends custom basis sets (from Basis Set Exchange) to ``self.footer``. Should be used in combination with the ``gen`` keyword.

        Args:
            name (str): name of basis set (look it up on Basis Set Exchange)
            add_all_elements (bool): whether the complete basis set should be added or just the elements of interest
            return_string (bool): if the basis set should be appended to the footer or returned as a string (no change to ``self``)

        Returns:
            nothing (if return_string is ``False``)
            string of basis set definition (if return string is ``True``)
        """
        import basis_set_exchange as bse
        assert isinstance(name, str), "need basis set name to be a string, for starters"

        try:
            basis_definition = ""
            if add_all_elements:
                basis_definition = bse.get_basis(name, fmt="gaussian94", header=False)
            else:
                elements = list(np.unique(self.get_molecule().atomic_numbers.view(np.ndarray)))
                basis_definition = bse.get_basis(name, fmt="gaussian94", header=False, elements=elements)

            if self.footer is None:
                self.footer = basis_definition
            else:
                self.footer += basis_definition
            self.footer += "\n"

        except Exception as e:
            raise ValueError(f"adding basis set {name} from basis set exchange failed!\n{e}")

    @classmethod
    def read_file(cls, filename, return_lines=False, extended_opt_info=False, fail_silently=True):
#    def read_fast(cls, filename, return_lines=False, extended_opt_info=False):
        """
        Reads a Gaussian``.out`` or ``.gjf`` file and populates the attributes accordingly.
        Only footers from ``opt=modredundant`` can be read automatically --  ``genecep`` custom basis sets, &c must be specified manually.

        Note:

        Will throw ``ValueError`` if there have been no successful iterations.

        Args:
            filename (str): path to the out file
            return_lines (Bool): whether the lines of the file should be returned
            extended_opt_info (Bool): if full parameters about each opt step should be collected
                (by default, only ``rms_displacement`` and ``rms_force`` are collected)
            fail_silently (Bool): if true, files that fail validation will just be omitted and parsing will continue.
                useful for monitoring jobs which are in-progress and may not have all properties written.
        Returns:
            ``GaussianFile`` object (or list of ``GaussianFile`` objects for Link1 files)
            (optional) the lines of the file (or list of lines of file for Link1 files) as Lines object
        """
        if re.search("gjf$", filename) or re.search("com$", filename):
            return cls._read_gjf_file(filename, return_lines)

        link1_lines = parse.split_link1_to_text(filename)
        files = []

        for link1idx, lines in enumerate(link1_lines):
            current_file = parse.read_file_fast(lines, filename, link1idx, extended_opt_info=extended_opt_info, fail_silently=fail_silently)
            if current_file is not None:
                files.append(current_file)

        if return_lines:
            link1_lines = parse.split_link1(filename)
            if len(link1_lines) == 1:
                return files[0], link1_lines[0]
            else:
                return files, link1_lines
        else:
            if len(link1_lines) == 1:
                return files[0]
            else:
                return files


