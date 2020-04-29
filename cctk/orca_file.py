import sys, re
import numpy as np

from abc import abstractmethod

from cctk import File, Molecule, ConformationalEnsemble
from cctk.helper_functions import get_symbol, compute_distance_between, compute_angle_between, compute_dihedral_between, get_number

import cctk.parse_orca as parse

class OrcaFile(File):
    """
    Generic class for all Orca `.inp` and `.out` files.

    Attributes:
        title (str): the title from the file
        ensemble (ConformationalEnsemble): `ConformationalEnsemble` instance
        header (str): keyword line or lines
        variables (dict): list of variables to specify (e.g. ``{"maxcore": 2000}``).
        blocks (dict): list of blocks to change specific settings
            In general, the key should be the block name and the value should be a list of desired lines.
            For instance, configuring a time-dependent DFT job might look like ``{"tddft": ["maxdim 5", "nroots 50"]}``
    """

    def __init__(self, ensemble, title=None, header=None, variables=None, blocks=None):
        if ensemble and isinstance(ensemble, ConformationalEnsemble):
            self.ensemble = ensemble
        if title and isinstance(title, str):
            self.title = title
        if header and isinstance(header, str):
            self.header = header

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
            energies = parse.read_energies()

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
    def write_molecule_to_file(cls, filename, molecule, header, variables=None, blocks=None, print_symbol=False):
        """
        Write a ``.inp``file using the given molecule.

        Args:
            filename (str): path to the new file
            molecule (Molecule): which molecule to use -- a``Molecule`` object.
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
        TODO: check indexing here

        Returns the last molecule (from an optimization job or other multi-molecule jobs) or the only molecule (from other jobs).

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

        return self.ensemble.molecules[num]
