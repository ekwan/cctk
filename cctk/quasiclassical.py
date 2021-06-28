"""
Functions to assist in sampling thermally excited states through quasiclassical approximations.
"""

import numpy as np
import math, copy, random

import cctk

"""
Constants:
"""

AMU_A2_FS2_PER_KCAL_MOL = 0.0004184

def get_quasiclassical_perturbation(molecule, temperature=298, return_velocities=False, which="quasiclassical"):
    """
    Perturbs molecule by treating each mode as a quantum harmonic oscillator and sampling from the distribution appropriate to the temperature.

    This is probably the only useful function in this file.

    Args:
        molecule (cctk.Molecule): molecule with vibrational modes
        temperature (float): temperature
        return velocities (bool): whether or not to return velocities
        which (str): ``classical`` or ``quasiclassical``

    Returns:
        new ``cctk.Molecule`` object
        energy above ground state (kcal/mol)
        velocities (cctk.OneIndexedArray)
    """
    assert isinstance(molecule, cctk.Molecule), "need a valid molecule"
    assert len(molecule.vibrational_modes) > 0, "molecule needs to have vibrational modes (try running a ``freq`` job)"

    assert isinstance(temperature, (int, float)), "temperature must be numeric"

    mol = copy.deepcopy(molecule)
    total_PE = 0
    total = 0

    velocities = np.zeros_like(molecule.geometry.view(np.ndarray)).view(cctk.OneIndexedArray)

    all_text = ""
    for mode in mol.vibrational_modes:
        PE, KE, TE, text = apply_vibration(mol, mode, temperature=temperature, which=which)
        total_PE += PE
        total += TE
        all_text += f"{text}\n"
        #### below code initializes velocities. pulling from jprogdyn Initializer lines 440 & onwards.

        if return_velocities:
            # mode velocity = sqrt(2 * KE / reduced mass) - want value in Å/fs. we randomize the sign
            # https://stackoverflow.com/questions/46820182/randomly-generate-1-or-1-positive-or-negative-integer
            mode_velocity = math.sqrt(2*KE*AMU_A2_FS2_PER_KCAL_MOL/mode.reduced_mass)
            mode_velocity *= (1 if random.random() < 0.5 else -1)

            for idx in range(1,molecule.num_atoms()+1):
                velocities[idx] += mode_velocity * mode.displacements[idx]

    if return_velocities:
        return mol, total_PE, total, text, velocities
    else:
        return mol, total_PE, total, text

def apply_vibration(molecule, mode, min_freq=50, temperature=298, verbose=False, which="quasiclassical"):
    """
    Apply a vibration to molecule ``molecule`` (modified in-place).

    Args:
        molecule (cctk.Molecule)
        mode (cctk.VibrationalMode)
        min_freq (float)
        temperature (float)
        verbose (bool)
        which (str): ``quasiclassical`` or ``classical``

    Returns:
        potential energy
        kinetic energy
        energy
        text
    """

    level = mode.choose_level(temperature)
    energy = mode.energy(level)

    shift = mode.random_displacement(level, method=which)
    max_shift = mode.classical_turning_point(level)
    rel_shift = shift/max_shift

    if mode.frequency < min_freq:
        rel_shift = 0

    molecule.geometry += mode.displacements * rel_shift * max_shift

    potential_energy = 0.5 * mode.force_constant * shift ** 2
    kinetic_energy = energy - potential_energy

    text = f"Mode {mode.frequency:.2f} ({mode.energy():.2f} kcal/mol)\t QC Level {level}\t Shift {rel_shift:.2%} of a potential {max_shift:.2f} Å\tPE = {potential_energy:.2f} kcal/mol\tk = {mode.force_constant:.2f} kcal/mol Å^-2"

    if verbose:
        print(text)

    return potential_energy, kinetic_energy, energy, text

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


