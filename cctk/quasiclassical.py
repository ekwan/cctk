"""
Functions to assist in sampling thermally excited states through quasiclassical approximations.
"""

import numpy as np
import math, copy

import cctk

def get_quasiclassical_perturbation(molecule, temperature=298):
    """
    Perturbs molecule by treating each mode as a quantum harmonic oscillator and sampling from the distribution appropriate to the temperature.

    This is probably the only useful function in this file.

    Args:
        molecule (cctk.Molecule): molecule with vibrational modes
        temperature (float): temperature

    Returns:
        new ``cctk.Molecule`` object
        energy above ground state (kcal/mol)
    """
    assert isinstance(molecule, cctk.Molecule), "need a valid molecule"
    assert len(molecule.vibrational_modes) > 0, "molecule needs to have vibrational modes (try running a ``freq`` job)"

    assert isinstance(temperature, (int, float)), "temperature must be numeric"

    mol = copy.deepcopy(molecule)
    total_PE = 0

    for mode in mol.vibrational_modes:
        PE, KE = apply_vibration(mol, mode, temperature=temperature)
        total_PE += PE

    return mol, total_PE

def apply_vibration(molecule, mode, min_freq=50, temperature=298, verbose=False):
    """
    Apply a vibration to molecule ``molecule`` (modified in-place).

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
    if verbose:
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


