import cctk
from cctk.helper_functions import get_symbol


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

    def write_file(self, filename, write_xyz=False, write_dir=None):
        """
        Write an SI file.

        Args:
            filename (str): path to the new file
            write_xyz (Bool): whether or not to write ``.xyz`` files for each molecule
            write_dir (str): where to write them too
        """
        first = True
        for title, (molecule, properties) in zip(self.titles, self.ensemble.items()):
            assert isinstance(molecule, cctk.Molecule), "molecule is not a valid Molecule object!"

            text = f"{title}\n"
            for key, value in generate_info(molecule, properties).items():
                text += f"{key}:\t{value}\n"

            text += f"Cartesian Coordinates (Å):\n"
            for index, Z in enumerate(molecule.atomic_numbers, start=1):
                line = molecule.get_vector(index)
                text += f"{get_symbol(Z):>2}       {line[0]:>13.6f} {line[1]:>13.6f} {line[2]:>13.8f}\n"

            text += "\n"

            if write_xyz and write_dir is not None:
                cctk.XYZFile.write_molecule_to_file(f"{write_dir}/{title}.xyz", molecule)

            if first:
                super().write_file(filename, text)
                first = False
            else:
                super().append_to_file(filename, text)


def generate_info(molecule, properties):
    info = {
        "Number of Atoms": molecule.num_atoms(),
        "Stoichiometry": molecule.formula(),
        "Charge": molecule.charge,
        "Multiplicity": molecule.multiplicity,
    }

    # for now manually handling route card and imaginaries, which typically aren't linked to cctk.Molecule.
    # long-term would be good to manually pass an extra info_dict from the calling environment
    # to avoid these ad hoc carveouts. ccw 3.8.21

    if "route_card" in properties:
        info["Route Card"] = properties["route_card"]

    if "imaginaries" in properties:
        info["Imaginary Frequencies (cm-1)"] = properties["imaginaries"]
    else:
        info["Imaginary Frequencies (cm-1)"] = "None"

    if "energy" in properties:
        info["Energy"] = properties["energy"]
    if "enthalpy" in properties:
        info["Enthalpy"] = properties["enthalpy"]
    if "gibbs_free_energy" in properties:
        info["Gibbs Free Energy"] = properties["gibbs_free_energy"]
    if "quasiharmonic_gibbs_free_energy" in properties:
        info["Gibbs Free Energy (Quasiharmonic Correction)"] = properties["quasiharmonic_gibbs_free_energy"]
    if "dipole_moment" in properties:
        info["Dipole Moment (Debye)"] = properties["dipole_moment"]

    return info

