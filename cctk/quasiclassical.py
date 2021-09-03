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
BOLTZMANN_CONSTANT = 0.001985875 # kcal/mol•Kn

def get_quasiclassical_perturbation(molecule, temperature=298, return_velocities=False, which="quasiclassical", mode_options=None):
    """
    Perturbs molecule by treating each mode as a quantum harmonic oscillator and sampling from the distribution appropriate to the temperature.

    This is probably the only useful function in this file.

    Args:
        molecule (cctk.Molecule): molecule with vibrational modes
        temperature (float): temperature
        return velocities (bool): whether or not to return velocities
        which (str): ``classical`` or ``quasiclassical``
        mode_options (dict):
            Options for how to initialize specific modes.
                key (int): 1-indexed number of vibrational mode (from smallest frequency to largest)
                val (dict):
                    velocity (str): one of "positive", "negative", "random", "zero"
                    displacement (bool): whether or not to displace

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

    if mode_options is None:
        mode_options = dict()

    all_text = ""
    for idx, mode in enumerate(mol.vibrational_modes):
        # enumerate is 0-indexed but GaussView, etc 1-index the modes. so we're 1-indexing here too.
        if idx+1 in mode_options:
            PE, KE, TE, mode_velocity, text = apply_vibration(mol, mode, temperature=temperature, which=which, **mode_options[idx+1])
        else:
            PE, KE, TE, mode_velocity, text = apply_vibration(mol, mode, temperature=temperature, which=which)
        total_PE += PE
        total += TE
        all_text += f"Mode {idx+1}: {text}\n"

        for idx in range(1,molecule.num_atoms()+1):
            velocities[idx] += mode_velocity * mode.displacements[idx]

    if return_velocities:
        return mol, total_PE, total, all_text, velocities
    else: # backwards compatibility
        return mol, total_PE, total, all_text

def apply_vibration(molecule, mode, min_freq=50, temperature=298, verbose=False, which="quasiclassical", displacement=True, velocity="random", **kwargs):
    """
    Apply a vibration to molecule ``molecule`` (modified in-place).

    Args:
        molecule (cctk.Molecule)
        mode (cctk.VibrationalMode)
        min_freq (float)
        temperature (float)
        verbose (bool)
        which (str): ``quasiclassical`` or ``classical``
        displacement (bool): whether or not to displace the mode
        velocity (str): ``positive``, ``negative``, ``random``, or ``zero``

    Returns:
        potential energy
        kinetic energy
        energy
        velocities
        text
    """

    if mode.frequency < 0:
        which = "ts"

    if which == "quasiclassical":
        level = mode.choose_level(temperature)
        energy = mode.energy(level)
        shift = mode.random_displacement(level=level, method=which)
        method = f"qc level {level}"
    elif which == "classical":
        energy = random_boltzmann_energy(temperature)
        shift = mode.random_displacement(energy=energy, method=which)
        method = "classical"
    elif which == "ts":
        energy = random_boltzmann_energy(temperature)
        shift = 0
        method = "ts"
    else:
        raise ValueError(f"``which`` must be ``classical``, ``quasiclassical``, or ``ts`` - {which} does not match!")

    # the rest is common to all methods

    # transition states and low-frequency modes do not get a starting displacement
    if not displacement or mode.frequency < min_freq:
        shift = 0

    max_shift = mode.classical_turning_point(energy=energy)
    if max_shift == 0.0:
        rel_shift = 0.0
        print("Warning: attempted to calculate relative shift when max shift is 0!")
    else:
        if shift > max_shift:
            print("Warning: requested shift of {shift:.4E} exceeds the max_shift of {max_shift:.4E}!")
            shift = max_shift
        rel_shift = shift/max_shift

    # apply displacements and compute energy breakdown
    molecule.geometry += mode.displacements * rel_shift * max_shift
    potential_energy = 0.5 * mode.force_constant * shift ** 2
    kinetic_energy = energy - potential_energy

    # mode velocity = sqrt(2 * KE / reduced mass) - want value in Å/fs
    # https://stackoverflow.com/questions/46820182/randomly-generate-1-or-1-positive-or-negative-integer
    mode_velocity = math.sqrt(2*kinetic_energy*AMU_A2_FS2_PER_KCAL_MOL/mode.reduced_mass)

    # choose velocity sign
    if velocity == "random":
        mode_velocity *= (1 if random.random() < 0.5 else -1)
    elif velocity == "negative":
        mode_velocity *= -1
    elif velocity == "zero":
        mode_velocity = 0
    elif velocity != "positive":
        raise ValueError(f"unknown value {velocity} for keywork ``velocity`` - must be ``positive``, ``negative``, ``random``, or ``zero``")

    text = f"{mode.frequency:.1f} cm-1 ({energy:4.2f} kcal/mol)\t{method}\t Shift {shift:5.2f} of {max_shift:4.2f} Å ({rel_shift:5.0%})"
    text += f"\tPE = {potential_energy:4.2f} kcal/mol\tKE = {kinetic_energy:4.2f} kcal/mol\tk = {mode.force_constant:.1f} kcal/mol Å^-2"
    if not displacement:
        text += "\n\t\tDisplacement manually set to zero!\n"
    if velocity == "zero":
        text += "\n\t\tVelocity manually set to zero!\n"
    if verbose:
        print(text)

    return potential_energy, kinetic_energy, energy, mode_velocity, text

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

def random_boltzmann_energy(temperature, cutoff=10, step1=0.01, step2=0.0001):
    """
    Randomly samples from the Boltzmann distribution appropriate for the given temperature.

    Arguments:
        temperature (int or float): in K
        cutoff (int or float): max energy considered, in kT
        step1: coarse numerical step, in kT
        step2: fine numerical step, in kT
    """
    kT = temperature * BOLTZMANN_CONSTANT

    random = np.random.uniform()

    # cumulative Boltzmann
    cumulative_boltzmann = lambda e: math.erf(math.sqrt(e))

    # now we need to numerically invert the cumulative Boltzmann
    # kT = 1 for all this math, we'll fix it at the end
    trial_energy = -1

    # scan up to cutoff kT, which should be more than enough
    for i in np.arange(0, cutoff, step1):
        if cumulative_boltzmann(i) > random:
            trial_energy = i - step1
            break

    if trial_energy == -1:
        return cutoff * kT

    # retry in smaller increments
    for i in np.arange(trial_energy, trial_energy+step1, step2):
        if cumulative_boltzmann(i) > random:
            trial_energy = i - step2
            break

    return trial_energy * kT
