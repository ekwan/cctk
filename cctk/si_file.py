import numpy as np

import cctk
from cctk.helper_functions import get_symbol, get_number


class SIFile(cctk.File):
    """
    Class representing Supporting Information files.

    Attributes:
        titles (list of str): title of each molecule
        ensemble (cctk.Ensemble): ``cctk.Ensemble`` of molecules to print
    """

    def __init__(self, ensemble, titles):
        if ensemble and isinstance(ensemble, cctk.Ensemble):
            self.ensemble = ensemble
        else:
            raise ValueError(f"invalid ensemble {ensemble}!")

        assert len(titles) == len(ensemble)
        self.titles = titles

    def write_file(self, filename, append=False):
        """
        Write an SI file.

        Args:
            filename (str): path to the new file
            append (Bool): whether or not to append to file
        """
        for title, (molecule, properties) in zip(self.titles, self.ensemble.items()):
            assert isinstance(molecule, cctk.Molecule), "molecule is not a valid Molecule object!"

            text = f"{title}\n"
            for key, value in generate_info(molecule, properties).items():
                text += f"{key}:\t{value}\n"

            for index, Z in enumerate(molecule.atomic_numbers, start=1):
                line = molecule.get_vector(index)
                text += f"{get_symbol(Z):>2}       {line[0]:>13.8f} {line[1]:>13.8f} {line[2]:>13.8f}\n"

            if append:
                text += "\n"
                super().append_to_file(filename, text)
            else:
                super().write_file(filename, text)


def generate_info(molecule, properties):
    info = {
        "Number of Atoms": molecule.num_atoms(),
        "Stoichiometry": molecule.formula(),
        "Charge": molecule.charge,
        "Multiplicity": molecule.multiplicity,
    }

    if "route_card" in properties:
        info["Route Card"] = properties["route_card"]

    if "energy" in properties:
        info["Energy"] = properties["energy"]
    if "enthalpy" in properties:
        info["Enthalpy"] = properties["enthalpy"]
    if "gibbs_free_energy" in properties:
        info["Gibbs Free Energy"] = properties["gibbs_free_energy"]
    if "quasiharmonic_gibbs_free_energy" in properties:
        info["Gibbs Free Energy (Quasiharmonic Correction)"] = properties["quasiharmonic_gibbs_free_energy"]

    return info

