import re
import numpy as np

from enum import Enum

from cctk import File, Molecule, ConformationalEnsemble
from cctk.helper_functions import get_symbol, get_corrected_free_energy

import cctk.parse_orca as parse

class OrcaJobType(Enum):
    """
    Class representing allowed Orca job types. Not an exhaustive list, but should be fairly comprehensive.

    The value should be the Orca keyword, to permit automatic assignment.

    All jobs have type ``SP`` by default.
    """

    SP = ("sp", ["energy", "scf_iterations"])
    """
    Single point energy calculation.
    """

    OPT = ("opt", ["rms_gradient", "rms_step", "max_gradient", "max_step"])
    """
    Geometry optimization.
    """

    FREQ = ("freq", ["gibbs_free_energy", "enthalpy", "frequencies", "temperature"])
    """
    Hessian calculation.
    """

    NMR = ("nmr", ["isotropic_shielding"])
    """
    NMR shielding prediction.
    """

    def __init__(self, value, expected_properties):
        self._value_ = value
        self.expected_properties = expected_properties

class OrcaFile(File):
    """
    Generic class for all Orca `.inp` and `.out` files.

    Attributes:
        ensemble (ConformationalEnsemble): `ConformationalEnsemble` instance
        job_types (list): list of ``OrcaJobType`` instances
        header (str): keyword line or lines
        variables (dict): list of variables to specify (e.g. ``{"maxcore": 2000}``).
        blocks (dict): list of blocks to change specific settings
            In general, the key should be the block name and the value should be a list of desired lines.
            For instance, configuring a time-dependent DFT job might look like ``{"tddft": ["maxdim 5", "nroots 50"]}``
        successful_terminations (int): number of successful terminations
        elapsed_time (float): total time for job in seconds
    """

    def __init__(self, job_types, ensemble=None,  header=None, variables=None, blocks=None):
        if job_types is not None:
            if not all(isinstance(job, OrcaJobType) for job in job_types):
                raise TypeError(f"invalid job type list {job_types}")
            self.job_types = job_types
        else:
            raise ValueError("need job types for new Orca file")

        if ensemble and isinstance(ensemble, ConformationalEnsemble):
            self.ensemble = ensemble
        else:
            self.ensemble = ConformationalEnsemble()

        if header and isinstance(header, str):
            self.header = header
        else:
            self.header = None

        if blocks and isinstance(blocks, dict):
            for lines in list(blocks.values()):
                assert isinstance(lines, list)
            self.blocks = blocks
        else:
            self.blocks = {}

        if variables and isinstance(variables, dict):
            self.variables = variables
        else:
            self.variables = {}

    @classmethod
    def read_file(cls, filename):
        if re.search("inp$", filename):
            return cls._read_inp_file(filename)

        multiple_lines = parse.split_multiple_inputs(filename)
        files = []

        for lines in multiple_lines:
            input_lines = parse.extract_input_file(lines)
            header = parse.read_header(input_lines)
            job_types = cls._assign_job_types(header)
            variables, blocks = parse.read_blocks_and_variables(input_lines)

            successful_scf_convergence = 0
            successful_opt = 0
            successful_freq = 0
            successful_NMR_EPR = 0
            is_scan_job = False
            # add identifiers for successful termination of other job types

            elapsed_time = 0
            for line in lines:
                if line.startswith("FINAL SINGLE POINT ENERGY"):                                        #### SCF converged at least once
                    successful_scf_convergence += 1
                elif line.strip().startswith("***        THE OPTIMIZATION HAS CONVERGED     ***"):      #### geometry converged
                    successful_opt += 1
                elif line.startswith("VIBRATIONAL FREQUENCIES"):                                        #### a frequency job was completed
                    successful_freq += 1
                elif line.startswith("Maximum memory used throughout the entire EPRNMR-calculation:"):  #### an EPR NMR job was completed
                    successful_NMR_EPR += 1
                elif line.strip().startswith("*    Relaxed Surface Scan    *"):                         #### this is a scan job
                    is_scan_job = True
                elif line.startswith("Sum of individual times         ..."):                            #### the job was completed  (note the '...' is key)
                    fields = line.split()
                    assert len(fields) == 9 or len(fields) == 10, f"unexpected number of fields on elapsed time line:\n{line}"
                    elapsed_time = float(fields[5])

            # different than G16 "successful termination"
            success = 0
            if successful_opt > 0:
                success += 1
            if successful_freq > 0:
                success += 1
            if successful_NMR_EPR > 0:
                success += 1
            if successful_scf_convergence > 0:
                success += 1

            energies, iters = parse.read_energies(lines)
            if len(energies) == 0:
                return None

            atomic_numbers, geometries = parse.read_geometries(lines, num_to_find=len(energies))
            assert len(geometries) >= len(energies), "can't have an energy without a geometry (cf. pigeonhole principle)"

            # this approach does not work with the option miniprint
            charge = lines.find_parameter("Total Charge           Charge          ....", 5, 4)[0]
            multip = lines.find_parameter("Multiplicity           Mult            ....", 4, 3)[0]

            #### TODO
            # detect Mayer bond orders

            f = OrcaFile(job_types, header=header, variables=variables, blocks=blocks)
            f.elapsed_time = elapsed_time
            f.successful_terminations = success

            molecules = [None] * len(geometries)
            properties = [{} for _ in range(len(geometries))]
            for idx, geom in enumerate(geometries):
                molecules[idx] = Molecule(atomic_numbers, geom, charge=charge, multiplicity=multip, bonds=None)
                if idx < len(energies):
                    properties[idx]["energy"] = energies[idx]
                properties[idx]["filename"] = filename
                properties[idx]["iteration"] = idx
                properties[idx]["scf_iterations"] = iters[idx]

            if multip > 1:
                s2 = lines.find_parameter("Expectation value of", 6, 5)
                for idx, spin_contam in enumerate(s2):
                    properties[idx]["S**2"] = spin_contam

            if OrcaJobType.OPT in job_types:
                rms_grad, max_grad, rms_step, max_step = parse.read_gradients(lines, len(properties))
                for idx in range(len(rms_grad)):
                    if idx < len(rms_grad):
                        properties[idx]["rms_gradient"] = rms_grad[idx]

                    if idx < len(max_grad):
                        properties[idx]["max_gradient"] = max_grad[idx]

                    if idx < len(rms_step):
                        properties[idx]["rms_step"] = rms_step[idx]

                    if idx < len(max_step):
                        properties[idx]["max_step"] = max_step[idx]

            if OrcaJobType.FREQ in job_types:
                properties[-1]["frequencies"] = sorted(parse.read_freqs(lines, successful_freq))

                enthalpies = lines.find_parameter("Total Enthalpy", expected_length=5, which_field=3)
                try:
                    properties[-1]["enthalpy"] = enthalpies[-1]
                except Exception:
                    pass

                gibbs = lines.find_parameter("Final Gibbs free", expected_length=7, which_field=5)
                try:
                    properties[-1]["gibbs_free_energy"] = gibbs[-1]
                except Exception:
                    pass

                try:
                    temperature = lines.find_parameter("Temperature", expected_length=4, which_field=2)
                    if len(temperature) > 0 and len(gibbs) > 0:
                        properties[-1]["temperature"] = temperature[-1]
                        corrected_free_energy = get_corrected_free_energy(gibbs[-1], properties[-1]["frequencies"],
                                                                      frequency_cutoff=100.0, temperature=temperature[-1])
                        properties[-1]["quasiharmonic_gibbs_free_energy"] = float(corrected_free_energy)
                except Exception:
                    pass

            if OrcaJobType.NMR in job_types:
                nmr_shifts = parse.read_nmr_shifts(lines, molecules[0].num_atoms())
                if nmr_shifts is not None:
                    properties[-1]["isotropic_shielding"] = nmr_shifts

            try:
                charges = parse.read_mulliken_charges(lines, successful_opt, is_scan_job)
                assert len(charges) == len(atomic_numbers)
                properties[-1]["mulliken_charges"] = charges
            except Exception:
                pass

            try:
                charges = parse.read_loewdin_charges(lines, successful_opt, is_scan_job)
                assert len(charges) == len(atomic_numbers)
                properties[-1]["lowdin_charges"] = charges
            except Exception:
                pass

            try:
                dipole = lines.find_parameter("Magnitude \(Debye\)", 4, 3)
                properties[-1]["dipole_moment"] = dipole[0]
            except Exception:
                pass

            for mol, prop in zip(molecules, properties):
                f.ensemble.add_molecule(mol, properties=prop)

            f.check_has_properties()
            files.append(f)

        if len(files) == 1:
            return files[0]
        else:
            return files

    @classmethod
    def _read_inp_file(cls, filename):
        print("reading ``.inp`` files is not currently supported :(")
        return None

    def write_file(self, filename, molecule=None, header=None, variables=None, blocks=None):
        """
        Write a ``.inp`` file, using object attributes. If no header is specified, the object's header will be used.

        Args:
            filename (str): path to the new file
            molecule (int): which molecule to use -- passed to ``self.get_molecule()``.
                Default is -1 (e.g. the last molecule), but positive integers will select from self.ensemble.molecules (0-indexed).
                A ``Molecule`` object can also be passed, in which case that molecule will be written to the file.
            header (str): header for new file
        """
        if molecule is None:
            molecule = -1
        if not isinstance(molecule, Molecule):
            molecule = self.ensemble.molecules[molecule]

        if header is None:
            header = self.header

        if variables is None:
            variables = self.variables

        if blocks is None:
            blocks = self.blocks

        self.write_molecule_to_file(filename, molecule, header, variables, blocks)

    @classmethod
    def write_molecule_to_file(cls, filename, molecule, header, variables={"maxcore": 2000}, blocks={"pal": ["nproc 16"]}, print_symbol=False):
        """
        Write an ``.inp`` file using the given molecule.

        Args:
            filename (str): path to the new file
            molecule (Molecule): which molecule to use -- a ``Molecule`` object.
            header (str): header for new file
            print_symbol (Bool): if atomic symbols should be printed instead of atomic numbers
        """
        assert isinstance(molecule, Molecule), "need a valid molecule to write a file!"
        assert isinstance(header, str), "can't write a file without a header"

        text = f"{header.strip()}\n"

        if variables is not None:
            assert isinstance(variables, dict), "blocks must be a dictionary"
            for k, v in variables.items():
                text += f"%{k} {v}\n"

        if blocks is not None:
            assert isinstance(blocks, dict), "blocks must be a dictionary"
            for k, v in blocks.items():
                text += f"%{k}\n"
                for line in v:
                    text += f"\t{line}\n"
                text += "end\n"

        text +="\n"

        text += f"* xyz {int(molecule.charge)} {int(molecule.multiplicity)}\n"
        for index, Z in enumerate(molecule.atomic_numbers, start=1):
            line = molecule.get_vector(index)
            if print_symbol:
                Z = get_symbol(Z)
                text += f"{Z:>2}       {line[0]:>13.8f} {line[1]:>13.8f} {line[2]:>13.8f}\n"
            else:
                text += f"{Z:2d}       {line[0]:>13.8f} {line[1]:>13.8f} {line[2]:>13.8f}\n"

        text += "*\n"
        text += "\n"

        #### write the file
        super().write_file(filename, text)

    def get_molecule(self, num=None):
        """
        Returns the last molecule (from an optimization job or other multi-molecule jobs) or the only molecule (from other jobs).

        If ``num`` is specified, returns that job (1-indexed for positive numbers). So ``job.get_molecule(3)`` will return the 3rd element of ``job.molecules``, not the 4th.
        """
        # some methods pass num=None, which overrides setting the default above
        if num is None:
            num = -1

        if not isinstance(num, int):
            raise TypeError("num must be int")

        return self.ensemble.molecule_list()[num]

    def num_imaginaries(self):
        """
        Returns the number of imaginary frequencies.
        """
        return len(self.imaginaries())

    def imaginaries(self):
        """
        Returns the imaginary frequencies, rounded to the nearest integer.
        """
        if (OrcaJobType.FREQ in self.job_types) and (self.ensemble[-1:,"frequencies"] is not None):
            freqs = self.ensemble[-1:,"frequencies"]
            if not isinstance(freqs, list) or len(freqs) == 0:
                return list()
            else:
                return list(map(int, np.array(freqs)[np.array(freqs) < 0]))
        else:
            return list()


    @classmethod
    def _assign_job_types(cls, header):
        """
        Assigns ``OrcaJobType`` objects from header. ``OrcaJobType.SP`` is assigned by default.

        Args:
            header (str): Orca header

        Returns:
            list of ``OrcaJobType`` objects
        """
        job_types = []
        for job_type in OrcaJobType:
            if re.search(f"{job_type.value}", str(header), re.IGNORECASE):
                job_types.append(job_type)
        if OrcaJobType.SP not in job_types:
            job_types.append(OrcaJobType.SP) # include SP in all jobs, whether specified in header or not
        if re.search("ScanTs", str(header), re.IGNORECASE):
            job_types.append(OrcaJobType.OPT) #ScanTs is an optimization job that does not contain the string "opt" in it's keyword
        return job_types

    def check_has_properties(self):
        """
        Checks that the file has all the appropriate properties for its job types, and raises ``ValueError`` if not.

        This only checks the last molecule in ``self.ensemble``, for now.
        """
        if self.successful_terminations > 0:
            if self.successful_terminations == 2 and ((OrcaJobType.OPT in self.job_types) and (OrcaJobType.FREQ in self.job_types)):
                return # opt and freq jobs should have three terminations
            for job_type in self.job_types:
                for prop in job_type.expected_properties:
                    if not self.ensemble.has_property(-1, prop):
                        raise ValueError(f"expected property {prop} for job type {job_type}, but it's not there!")
        else:
            return


