from cctk import XYZFile, ConformationalEnsemble, GaussianFile
import numpy as np

#### Usage: ``python generate_conformers.py``

output_file = XYZFile.read_file('Ac-F2Ala-F2Ala-OMe.xyz')
output_file.molecule.assign_connectivity()

angles = range(0, 360, 120)
bin_angles = range(0, 360, 180)

to_rotate = [[1, 3, 5, 7], [9, 11, 13, 15], [5, 3, 6, 8], [12, 11, 14, 16]]
to_bin_rotate = [[27, 26, 1, 2], [8, 6, 9, 10], [16, 14, 17, 18]]

#### now to employ some recursion...
def rotate_angles_one_by_one (idx, angles, thetas, structures):
    if idx >= len(angles):
        return structures

    else:
        new_structures = [None] * (len(structures) * len(thetas))
        current_idx = 0
        for structure in structures:
            for theta in thetas:
                new_structures[current_idx] = structure.set_dihedral(*angles[idx], theta, check_result=False)
                current_idx += 1

    print(len(structures))
    return rotate_angles_one_by_one(idx+1, angles, thetas, new_structures)

mols = rotate_angles_one_by_one(0, to_rotate, angles, [output_file.molecule])
mols = rotate_angles_one_by_one(0, to_bin_rotate, bin_angles, mols)

for idx, molecule in enumerate(mols):
    try:
        molecule.check_for_conflicts()
        GaussianFile.write_molecule_to_file(f"conformer_dft_{idx:05d}.gjf", molecule, "#p opt b3lyp/6-31g(d) empiricaldispersion=gd3bj scrf=(smd,solvent=diethylether)", None)
    except:
        pass
