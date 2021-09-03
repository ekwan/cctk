import numpy as np
import sys

sys.path.insert(0,'/Users/cwagen/code/cctk')

import cctk
import cctk.quasiclassical as qc

#### This is a script to let you specify vibrational excited states manually, to compare with values from other programs.
#### This was used to test against Jprogdyn.
#### Corin Wagen, July 2020

path = "test/static/methane_hpmodes.out"
file = cctk.GaussianFile.read_file(path)

mol = file.get_molecule()

total_PE = 0

for mode in mol.vibrational_modes:
    print()
    print(mode)

    level = int(input("level: "))
    energy = mode.energy(level)

    rel_shift = float(input("rel_shift (%): ")) / 100
    max_shift = mode.classical_turning_point(energy=energy)

    shift = rel_shift * max_shift

    mode_coords = mode.displacements
    mol.geometry += mode.displacements * rel_shift * max_shift

    potential_energy = 0.5 * mode.force_constant * shift ** 2
    total_PE += potential_energy

    print(f"Mode {mode.frequency:.2f} ({mode.energy():.2f} kcal/mol)\t QC Level {level}\t Shift {rel_shift:.2%} of a potential {max_shift:.2f} Å\tPE = {potential_energy:.2f} kcal/mol\tk = {mode.force_constant:.2f} kcal/mol Å^-2")

cctk.GaussianFile.write_molecule_to_file("test/static/methane_perturbed.gjf", mol, route_card="#p sp hf/sto-3g")
