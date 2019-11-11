import numpy as np

from cctk import XYZFile, OutputFile
from cctk.helper_functions import get_number


class XYZData (OutputFile):
    """
    Creates output file instances of the specific type through factory methods.
    """ 

    @classmethod
    def read_xyz(cls, filename):
        lines = super().read_file(filename)
        num_atoms = 0

        try:
            num_atoms = int(lines[0])
        except:
            raise ValueError("can't get the number of atoms from the first line!")

        assert num_atoms == (len(lines) - 2), "wrong number of atoms!"

        title = lines[1]

        atoms = [None] * num_atoms
        geometry = [None] * num_atoms

        for index, line in enumerate(lines[2:]):
            pieces = list(filter(None, line.split(" ")))
            try:
                atoms[index] = get_number(pieces[0])
                geometry[index] = [float(pieces[1]), float(pieces[2]), float(pieces[3])]
            except:
                raise ValueError(f"can't parse line {index+2}!")

        return XYZFile(atoms, geometry, title)

