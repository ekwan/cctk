import re
import numpy as np

import cctk
from cctk.helper_functions import get_symbol, get_number


class XYZFile(cctk.File):
    """
    Class representing plain ``.xyz`` files.

    Attributes:
        title (str): the title from the file
        molecule (Molecule): `Molecule` instance
    """

    def __init__(self, molecule, title=None):
        if molecule and isinstance(molecule, cctk.Molecule):
            self.molecule = molecule
        if title and (isinstance(title, str)):
            self.title = title

    @classmethod
    def read_file(cls, filename, charge=0, multiplicity=1):
        """
        Factory method to create new XYZFile instances.


        Arguments:
            filename (str): path to ``.xyz`` file
            charge (int): charge of resultant molecule
            multiplicity (int): multiplicity of resultant molecule
        """
        lines = super().read_file(filename)
        xyzfile = cls.file_from_lines(lines)

        assert isinstance(charge, int), "charge must be integer"
        assert isinstance(multiplicity, int), "multiplicity must be integer"
        assert multiplicity > 0, "multiplicity must be a positive integer"

        xyzfile.molecule.charge = charge
        xyzfile.molecule.multiplicity = multiplicity

        return xyzfile

    @classmethod
    def file_from_lines(cls, lines):
        num_atoms = 0
        try:
            num_atoms = int(lines[0])
        except:
            raise ValueError("can't get the number of atoms from the first line!")

        title = lines[1]

        atomic_numbers = np.zeros(shape=num_atoms, dtype=np.int8)
        geometry = np.zeros(shape=(num_atoms, 3))

        for index, line in enumerate(lines[2:]):
            # ignore blank lines
            if len(line.strip()) == 0:
                continue

            pieces = list(filter(None, line.split(" ")))
            try:
                if re.match("[0-9]", pieces[0]):
                    atomic_numbers[index] = int(pieces[0])
                elif re.match("([A-Za-z])+([0-9])+", pieces[0]):
                    # mdtraj writes in this format, for some reason
                    m = re.match("([A-Za-z])+([0-9])+", pieces[0])
                    atomic_numbers[index] = int(get_number(m.group(1)))
                else:
                    atomic_numbers[index] = int(get_number(pieces[0]))
                geometry[index][0] = float(pieces[1])
                geometry[index][1] = float(pieces[2])
                geometry[index][2] = float(pieces[3])
            except:
                raise ValueError(f"can't parse line {index+2}: {line}")

        assert num_atoms == len(atomic_numbers), "wrong number of atoms!"
        molecule = cctk.Molecule(atomic_numbers, geometry)
        return XYZFile(molecule, title)

    @classmethod
    def write_molecule_to_file(cls, filename, molecule, title="title", append=False):
        """
        Write an ``.xyz`` file, using object attributes.

        Args:
            filename (str): path to the new file
            molecule (Molecule): molecule to write
            title (str): title of file
            append (Bool): whether or not to append to file
        """
        assert isinstance(molecule, cctk.Molecule), "molecule is not a valid Molecule object!"

        text = f"{molecule.num_atoms()}\n"
        text += f"{title}\n"

        for index, Z in enumerate(molecule.atomic_numbers, start=1):
            line = molecule.get_vector(index)
            text += f"{get_symbol(Z):>2}       {line[0]:>13.8f} {line[1]:>13.8f} {line[2]:>13.8f}\n"

        if append:
            super().append_to_file(filename, text)
        else:
            super().write_file(filename, text)

    def write_file(self, filename):
        """
        Write an ``.xyz`` file, using object attributes.
        """
        self.write_molecule_to_file(filename, self.molecule, title=self.title)

    @classmethod
    def read_trajectory(cls, filename):
        """
        Read an ``.xyz`` trajectory file, which is just a bunch of concatenated ``.xyz`` files.
        Currently the files must be separated by nothing (no blank line, just one after the other) although this may be changed in future.

        Args:
            filename (str): path to file

        Returns:
            list of ``cctk.XYZFile`` objects in the order they appear in the file
        """
        files = []
        lines = super().read_file(filename)

        current_lines = list()
        for line in lines:
            if re.search(r"^\s*\d+$", line):
                if len(current_lines) > 0:
                    files.append(cls.file_from_lines(current_lines))
                    current_lines = list()
            current_lines.append(line)

        return files

    @classmethod
    def read_ensemble(cls, filename, conformational=False):
        """
        Alias for read_trajectory.
        """
        files = cls.read_trajectory(filename)

        ensemble = None
        if conformational:
            ensemble = cctk.ConformationalEnsemble()
        else:
            ensemble = cctk.Ensemble()

        for f in files:
            ensemble.add_molecule(f.molecule)

        return ensemble

    @classmethod
    def write_ensemble_to_file(cls, filename, ensemble, title=None):
        """
        Write a ``cctk.Ensemble`` to a single ``.xyz`` file. Can be viewed in MOLDEN.
        """
        assert isinstance(ensemble, cctk.Ensemble), f"ensemble {ensemble} is not a cctk.Ensemble"

        if title is None:
            title = "title"
        if isinstance(title, str):
            title = [title for _ in range(len(ensemble))]
        assert len(title) == len(ensemble)

        for idx, (molecule, title) in enumerate(zip(ensemble._items, title)):
            if idx == 0:
                cls.write_molecule_to_file(filename, molecule, title=title, append=False)
            else:
                cls.write_molecule_to_file(filename, molecule, title=title, append=True)
