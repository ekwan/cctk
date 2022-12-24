import re, warnings
import numpy as np

import cctk
from cctk.helper_functions import get_symbol, get_number


class XYZFile(cctk.File):
    """
    Class representing plain ``.xyz`` files.

    Attributes:
        titles (list of str): the title or titles from the file
        ensemble (Ensemble): `Ensemble` instance
        molecule (Molecule): `Molecule` instance representing the first molecule in the file. deprecated, but present for backwards compatibility.
    """

    def __init__(self, ensemble, titles):
        assert isinstance(ensemble, cctk.Ensemble), "ensemble must be cctk.Ensemble"
        self.ensemble = ensemble

        # backwards compatibility
        self.molecule = ensemble.molecule_list()[0]

        assert isinstance(titles, list), "title must be list"
        self.titles = titles

    def __getattribute__(self, name):
        if name == "molecule":
            warnings.warn("XYZFile attribute ``molecule`` will be removed in upcoming releases of cctk. Use ``ensemble`` attribute instead!", DeprecationWarning, stacklevel=2)
        return cctk.File.__getattribute__(self, name)

    @classmethod
    def read_file(cls, filename, charge=0, multiplicity=1, conformational=False):
        """
        Factory method to create new XYZFile instances.

        Arguments:
            filename (str): path to ``.xyz`` file
            charge (int): charge of resultant molecule
            multiplicity (int): multiplicity of resultant molecule
            conformational (bool): whether or not it's a conformational ensemble
        """
        assert isinstance(charge, int), "charge must be integer"
        assert isinstance(multiplicity, int), "multiplicity must be integer"
        assert multiplicity > 0, "multiplicity must be a positive integer"

        ensemble = cctk.Ensemble()
        if conformational:
            ensemble = cctk.ConformationalEnsemble()
        titles = list()

        lines = super().read_file(filename)
        current_lines = list()
        for line in lines:
            if re.search(r"^\s*\d+$", line) and len(current_lines) > 2:
                if len(current_lines) > 0:
                    t, m = cls.mol_from_lines(current_lines, charge=charge, multiplicity=multiplicity)
                    ensemble.add_molecule(m)
                    titles.append(t)
                    current_lines = list()
            current_lines.append(line)

        # catch the last molecule
        if len(current_lines) > 0:
            t, m = cls.mol_from_lines(current_lines, charge=charge, multiplicity=multiplicity)
            ensemble.add_molecule(m)
            titles.append(t)

        return XYZFile(ensemble, titles)

    @classmethod
    def mol_from_lines(cls, lines, charge=0, multiplicity=1):
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
        molecule = cctk.Molecule(atomic_numbers, geometry, charge=charge, multiplicity=multiplicity)
        return title, molecule

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

    def write_file(self, filename, idx=-1):
        """
        Write an ``.xyz`` file, using object attributes.

        Args:
            idx (int): the index of the molecule to write
        """
        assert isinstance(idx, int), "idx must be int"
        self.write_molecule_to_file(filename, self.get_molecule(idx), title=self.titles[idx])

    @classmethod
    def read_trajectory(cls, filename, **kwargs):
        """
        Post refactoring, just an alias for ``XYZFile.read_file()``.
        """
        return cls.read_file(filename, **kwargs)

    @classmethod
    def read_ensemble(cls, filename, **kwargs):
        """
        Post refactoring, just an alias for ``XYZFile.read_file()``.
        """
        return cls.read_file(filename, **kwargs)

    @classmethod
    def write_ensemble_to_file(cls, filename, ensemble, titles=None):
        """
        Write a ``cctk.Ensemble`` to a single ``.xyz`` file. Can be viewed in MOLDEN.

        Arguments:
            filename (str): path to ``.xyz`` file
            ensemble (Ensemble): the collection of structures to write
            titles (None, str, or list of str): if None, the titles of the Molecules will be used;
                                                if one str, then the same name will be used for all
                                                molecules; if iterable, then the titles are assumed
                                                to parallel the indexing of the list
        """
        assert isinstance(filename, str), f"got {type(filename)} for filename but expected str"
        assert isinstance(ensemble, cctk.Ensemble), f"ensemble {ensemble} is not a cctk.Ensemble"
        assert len(ensemble)>0, "can't write empty Ensemble to xyz file"

        if titles is None:
            pass
        elif isinstance(titles, str):
            assert len(titles) > 0, "zero length title not allowed"
            titles = [titles] * len(ensemble)
        elif isinstance(title, (list,np.ndarray)):
            assert len(titles) == len(ensemble)
            for i,title in enumerate(titles):
                assert isinstance(title, str), f"got {type(filename)} at index {i} of titles, but expected str"
        else:
            raise ValueError(f"got {type(titles)} for title but expected None, str, or iterable")

        for idx,molecule in enumerate(ensemble._items):
            append = idx > 0
            if titles is None:
                title = molecule.name
                if not (isinstance(title, str) and len(title)>0):
                    title = "title"
            else:
                title = titles[idx]
            cls.write_molecule_to_file(filename, molecule, title=title, append=append)

    def get_molecule(self, num=None):
        """
        Returns a given molecule.

        If ``num`` is specified, returns ``self.ensemble.molecule_list()[num]``
        """
        # some methods pass num=None, which overrides setting the default above
        if num is None:
            num = -1
        assert isinstance(num, int), "num must be int"
        return self.ensemble.molecule_list()[num]


