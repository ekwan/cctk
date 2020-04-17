.. _tutorial_03:

=======================================
Tutorial 03: Bulk Analysis/Resubmission
=======================================

Objectives
==========

This tutorial will teach:

- Creating command-line *cctk* scripts.
- Reading/writing data from output files. 

Overview
========

*cctk* comes with several scripts which enable rapid analysis and resubmission of large numbers of Gaussian output files. 

In this tutorial, we will dissect these scripts and explain how exactly they work. 

Monitoring a Single Job
=======================

To monitor a single job, we want to be able to see the energy, forces, and SCF convergence for each iteration. 
We first read in the output file as a ``GaussianFile`` object and extract the relevant parameters::

    filename = sys.argv[1]
    print(f"reading file {filename}")
    try:
        output_file = GaussianFile.read_file(filename)

        energies = output_file.energies
        scf_iter = output_file.scf_iterations
        rms_displacements = output_file.rms_displacements
        rms_forces = output_file.rms_forces

If we are mid-optimization, we will usually have uneven numbers of energies, displacements, and forces. 
We can fix this problem by simply discarding the values from the partially completed iteration::

    if len(energies) > len(rms_forces):
        #### sometimes you just catch the job in the middle
        if len(rms_forces) == 0:
            raise ValueError("not all necessary parameters calculated yet")

        energies = energies[:len(rms_forces)-1]
        scf_iter = scf_iter[:len(rms_forces)-1]
        rms_displacements = rms_displacements[:len(rms_forces)-1]

We can also print extra information for a successful job::

    if output_file.succesful_terminations:
        print("Optimization converged!")
        print(f"{output_file.num_imaginaries()} imaginary frequencies")

We may also want to visualize how bond distances are changing. 
We don't know *a priori* which bond distances will be interesting, but choosing the largest atom and its largest neighbor seems reasonable. 
Here, we use the ``get_geometric_parameters`` function to output the distances from each iteration::

    mol = output_file.molecules[0]
    biggest_atom = np.argmax(mol.atomic_numbers) + 1
    nearby_atoms = mol.get_adjacent_atoms(biggest_atom)
    biggest_nearby_atom = nearby_atoms[np.argmax([mol.atomic_numbers[x] for x in nearby_atoms])]
        
    ... 

    distances = output_file.molecules.get_geometric_parameters('distance',biggest_atom,biggest_nearby_atom)

All that remains is to print the output in an organized fashion::

    print("{0:5} {1:20} {2:20} {3:15} {4:15} {5:20} {6:15}".format(
        "#",
        "Energy (Hartree)",
        "Rel Energy (kcal)",
        "SCF Cycles",
        "RMS Force",
        "RMS Displacement",
        f"Distance({mol.atom_string(biggest_atom)},{mol.atom_string(biggest_nearby_atom)})"
    ))

    ...

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

The full script (``monitor.py``) is shown below::

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

        if output_file.succesful_terminations:
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
