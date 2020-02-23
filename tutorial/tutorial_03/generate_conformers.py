import numpy as np
import copy

from cctk import XYZFile, ConformationalEnsemble, GaussianFile

#### Usage: ``python generate_conformers.py``
#### This script takes an input ``.xyz`` file and outputs ~600 different conformations. 
#### Corin Wagen and Eugene Kwan, 2019

output_file = XYZFile.read_file('Ac-F2Ala-F2Ala-OMe.xyz')
output_file.molecule.assign_connectivity()

#### here we define the different choices for each angle.
angles = range(0, 360, 120)
bin_angles = range(0, 360, 180)

#### here we define all the rotatable bonds, and the bonds with two conformers (like amides)
to_rotate = [[1, 3, 5, 7], [9, 11, 13, 15], [5, 3, 6, 8], [12, 11, 14, 16]]
to_bin_rotate = [[27, 26, 1, 2], [8, 6, 9, 10], [16, 14, 17, 18]]

#### now to employ some recursion...
def rotate_angles_one_by_one (idx, angles, thetas, structures):
    """
    This script takes a set of structures, and outputs a (longer) set of structures where the given bond has been rotated.

    Args:
        idx (int): the current position in ``angles``
        angles (list of 4-element lists): the list of angles to recurse through and adjust
        thetas (list of float): the list of dihedral angles to set each bond to
        structures (list of cctk.Molecule): the current list of structures

    Returns:
        list of cctk.Molecule objects (with len(thetas) * len(structures) elements)
    """
    if idx >= len(angles):
        return structures

    else:
        new_structures = [None] * (len(structures) * len(thetas))
        current_idx = 0
        for structure in structures:
            for theta in thetas:
                new_structures[current_idx] = copy.deepcopy(structure.set_dihedral(*angles[idx], theta, check_result=False))
                current_idx += 1

    print(len(structures))
    return rotate_angles_one_by_one(idx+1, angles, thetas, new_structures)

mols = rotate_angles_one_by_one(0, to_rotate, angles, [output_file.molecule])
mols = rotate_angles_one_by_one(0, to_bin_rotate, bin_angles, mols)

for idx, molecule in enumerate(mols):
    try:
        molecule.check_for_conflicts()
        GaussianFile.write_molecule_to_file(f"conformer_{idx:05d}.gjf", molecule, "#p opt pm7", None)
    except:
        pass
