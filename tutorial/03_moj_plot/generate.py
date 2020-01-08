import numpy as np
from cctk import ConformationalEnsemble, GaussianFile

path = "starting_structure.gjf"
cent = 1
lg = 7
nu = 8

file = GaussianFile.read_file(path)
mol = file.get_molecule()
mol.assign_connectivity()

footer = f"B {cent} {lg} F\nB {cent} {nu} F\n"

for lg_dist in np.arange(1.5, 3.3, 0.1):
    for nu_dist in np.arange(1.2, 3.0, 0.1):
        mol.set_distance(cent, lg, lg_dist)
        mol.set_distance(cent, nu, nu_dist)

        mol.check_for_conflicts()
        GaussianFile.write_molecule_to_file(
            f"scan_{int(round(lg_dist*100))}_{int(round(nu_dist*100))}.gjf",
            mol,
            "#p opt=modredundant b3lyp/6-31+g(d) scrf=(smd, solvent=tetrahydrofuran)",
            footer=footer
        )
