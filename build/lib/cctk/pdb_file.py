from cctk import File
from cctk.helper_functions import get_symbol

class PDBFile(File):
    """
    Generic class for all ``.pdb`` files.
    """

    def __init__(self, molecule, title=None):
        pass

    @classmethod
    def read_file(cls, filename):
        pass

    @classmethod
    def write_molecule_to_file(cls, filename, molecule, num=1, append=False):
        """
        Write a ``.pdb`` file, using object attributes.

        Args:
            filename (str): path to the new file
            molecule (Molecule): ``Molecule`` object
            num (int): model number
            append (Bool): whether to write to file normally or append
        """
        text = f"MODEL {num}\n"

        for idx, Z in enumerate(molecule.atomic_numbers, start=1):
            line = molecule.get_vector(idx)
            symb = get_symbol(Z).upper()
            text += f"HETATM {idx:>4}  {symb:<2}    *     0     {line[0]:7.3f} {line[1]:7.3f} {line[2]:7.3f}  1.00  0.00          {symb:>2}\n"

        text += f"ENDMDL\n"

        if append:
            super().append_to_file(filename, text)
        else:
            super().write_file(filename, text)


    @classmethod
    def write_ensemble_to_trajectory(cls, filename, ensemble):
        """
        Writes a ``ConformationalEnsemble`` to a trajectory file.

        Args:
            filename (str): where to write the file
            ensemble (Ensemble): ``Ensemble`` object to write
        """
        for idx, molecule in enumerate(ensemble.molecules):
            if idx == 0:
                cls.write_molecule_to_file(filename, molecule, num=idx+1, append=False)
            else:
                cls.write_molecule_to_file(filename, molecule, num=idx+1, append=True)

