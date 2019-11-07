import sys
import os
import numpy as np

sys.path.append(os.path.relpath('../cctk'))

from cctk import GaussianData, Molecule, GaussianJob

output_file = GaussianData.read_opt('cctk/scripts/thiourea_catalyst.out')

if output_file.success == 0: 
    print("termination not successful!")

molecule = Molecule(output_file.atoms, output_file.get_final_geometry())
print(molecule.formula())
molecule.assign_connectivity()

if not os.path.exists('cctk/scripts/scan_angle'):
    os.makedirs('cctk/scripts/scan_angle')

angles = [0, 30, 60, 90, 120, 180]
molecule.set_dihedral(43, 27, 10, 12, 0)
for angle in angles:    
    molecule.set_dihedral(49, 9, 27, 43, angle)

    chlorines = molecule.get_atoms_by_symbol('Cl')
    for chlorine in chlorines: 
        molecule.remove_atom(chlorine)

    number = molecule.add_atom_at_centroid('Cl', [42, 43])
    molecule.set_distance(9, number, 2.5)

    input_file = GaussianJob.create_constrained_opt(molecule.atoms, molecule.geometry, [['D', 49, 9, 27, 43, 'F']])
    input_file.write_file(filename=f"cctk/scripts/scan_angle/catalyst_{angle}.gjf")
