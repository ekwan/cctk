"""
Functions to assist in sampling thermally excited states through quasiclassical approximations.
"""

import numpy as np
import math, copy

import cctk

def get_quasiclassical_vibration(molecule, temperature=298):
    assert isinstance(molecule, cctk.Molecule), "need a valid molecule"

    mol = copy.deepcopy(molecule)
    total_PE = 0

    for mode in mol.vibrational_modes:
        PE, KE = apply_vibration(mol, mode, temperature=temperature)
        total_PE += PE

    print(total_PE)

    return mol

def apply_vibration(molecule, mode, min_freq=50, temperature=298):
    """
    Args:
        molecule (cctk.Molecule)
        mode (cctk.VibrationalMode)

    Returns:
        potential energy
        kinetic energy
    """

    level = mode.choose_level(temperature)
    energy = mode.energy(level)

    shift = mode.random_displacement(level)
    max_shift = mode.classical_turning_point(level)
    rel_shift = shift/max_shift

    if mode.frequency < min_freq:
        rel_shift = 0

    mode_coords = mode.displacements
    molecule.geometry += mode.displacements * rel_shift * max_shift

    potential_energy = 0.5 * mode.force_constant * shift ** 2

    kinetic_energy = energy - potential_energy

    # here we could compute atom velocities if we wanted to! initializer lines 440-480

    print(f"Mode {mode.frequency:.2f} ({mode.energy():.2f} kcal/mol)\t QC Level {level}\t Shift {rel_shift:.2%} of a potential {max_shift:.2f} Å\tPE = {potential_energy:.2f} kcal/mol\tk = {mode.force_constant:.2f} kcal/mol Å^-2")

    return potential_energy, kinetic_energy


def get_hermite_polynomial(n):
    """
    Returns a ``np.poly1d`` object representing the degree-n Hermite polynomial.

    Adapted from https://scipython.com/blog/the-harmonic-oscillator-wavefunctions/.
    """
    assert isinstance(n, int) and n >= 0, "need positive integer"

    Hr = [None] * (n + 1)
    Hr[0] = np.poly1d([1.,])

    if n > 0:
        Hr[1] = np.poly1d([2., 0.])

    if n > 1:
        for v in range(2, n+1):
            Hr[v] = Hr[1]*Hr[v-1] - 2*(v-1)*Hr[v-2]
    return Hr[n]


