import sys
import numpy as np
from cctk import GaussianFile, Molecule

#### This is a script to monitor the output of Gaussian files. 
#### In contrast to ``analyze.py``, this script analyzes only a *single* file in depth. 
#### If the file has not successfully achieved SCF convergence at least once, this file will not display any information. 
#### usage: ``python monitor.py path/to/output.out``

#### Corin Wagen and Eugene Kwan, 2019

filename = sys.argv[1]
print(f"reading file {filename}")
try:
    output_file = GaussianFile.read_file(filename)

    energies = output_file.energies
    scf_iter = output_file.scf_iterations
    rms_displacements = output_file.rms_displacements
    rms_forces = output_file.rms_forces

    if len(energies) > len(rms_forces):
        #### sometimes you just catch the job in the middle
        if len(rms_forces) == 0:
            raise ValueError("not all necessary parameters calculated yet")

        energies = energies[:len(rms_forces)-1]
        scf_iter = scf_iter[:len(rms_forces)-1]
        rms_displacements = rms_displacements[:len(rms_forces)-1]

    if output_file.success:
        print("Optimization converged!")
        print(f"{output_file.num_imaginaries()} imaginary frequencies")

    #### often you care about the largest atom and its neighbors... so this will automatically print that bond distance 
    mol = output_file.molecules[0]
    biggest_atom = np.argmax(mol.atomic_numbers) + 1
    nearby_atoms = mol.get_adjacent_atoms(biggest_atom)
    biggest_nearby_atom = nearby_atoms[np.argmax([mol.atomic_numbers[x] for x in nearby_atoms])]

    print("{0:5} {1:20} {2:20} {3:15} {4:15} {5:20} {6:15}".format(
        "#",
        "Energy (Hartree)",
        "Rel Energy (kcal)",
        "SCF Cycles",
        "RMS Force",
        "RMS Displacement",
        f"Distance({mol.atom_string(biggest_atom)},{mol.atom_string(biggest_nearby_atom)})"
    ))

    distances = output_file.molecules.get_geometric_parameters('distance',biggest_atom,biggest_nearby_atom)

    min_energy = np.min(energies)
    for i, energy in enumerate(energies):
        rel_energy = (energy - min_energy) * 627.509 #### convert to kcal
        print("{0:5d} {1:20.5f} {2:20.5f} {3:15d} {4:15.5f} {5:20.5f} {6:15.3f}".format(
            i+1,
            energy,
            rel_energy,
            scf_iter[i],
            rms_forces[i],
            rms_displacements[i],
            distances[i]
        ))

except ValueError as e:
    print(f"job has not finished any iterations: {e}")
