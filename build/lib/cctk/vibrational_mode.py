import math
import numpy as np

import cctk
from cctk.quasiclassical import get_hermite_polynomial

# constants
MAX_QHO_LEVEL = 10000
MIN_FREQUENCY = 2
MIN_TEMPERATURE = 10
MAX_ZPE_RATIO = 0.999999

BOLTZMANN_CONSTANT = 0.001985875 # kcal/mol•K

class VibrationalMode:
    """
    Most code adapted from ``jprogdyn``. Displacements will be very low accuracy unless ``freq=hpmodes`` is enabled.

    Values from Gaussian, for now: see https://gaussian.com/vib/.

    Attributes:
        frequency (float): frequency, in cm-1
        force_constant (float): force constant, in kcal/mol per Å
        reduced_mass (float): mass, in amus
        intensity (float): IR intensity
        displacements (cctk.OneIndexedArray): atom displacements
        velocities (cctk.OneIndexedArray): atom velocities

    """
    def __init__(self, frequency, force_constant, reduced_mass, intensity, displacements):
        assert isinstance(frequency, float)
        self.frequency = frequency

        assert isinstance(force_constant, float)
        self.force_constant = force_constant

        assert isinstance(reduced_mass, float)
        self.reduced_mass = reduced_mass

        assert isinstance(intensity, float)
        self.intensity = intensity

        assert isinstance(displacements, cctk.OneIndexedArray)
        self.displacements = displacements

    def __str__(self):
        return f"Vibrational mode ({self.frequency:.2f} cm-1, {self.reduced_mass:.2f} amus, {self.force_constant:.2f} kcal/mol Å**-2)"

    def __repr__(self):
        return f"Vibrational mode ({self.frequency:.2f} cm-1, {self.reduced_mass:.2f} amus, {self.force_constant:.2f} kcal/mol Å**-2)"

    def choose_level(self, temperature=298):
        if temperature < MIN_TEMPERATURE:
            return 0

        # zpe_ratio is probability of being in level i vs level i+1, by quantum harmonic oscillator
        zpe_ratio = math.exp( -2 * self.energy() / (BOLTZMANN_CONSTANT * temperature))
        if zpe_ratio > MAX_ZPE_RATIO:
            zpe_ratio = MAX_ZPE_RATIO

        # probability of being in state 0 is equal to 1 - zpe_ratio
        # 1 = P(0) + P(1) + P(2) + ... = P + P * zpe_ratio + P * zpe_ratio ** 2 + ...
        # 1 = P(0) / (1 - zpe_ratio) bc geometric series
        P = 1.0 - zpe_ratio

        random = np.random.uniform()
        level = 0
        while level < MAX_QHO_LEVEL:
            if random < P:
                return level
            else:
                P += P * zpe_ratio
                level += 1

        return level

    def energy(self, level=0):
        """
        Calculate energy as a function of level. By default returns zero-point energy (level = 0).

        Args:
            level (int): which vibrational level the mode is in

        Returns:
            energy (kcal/mol)
        """
        assert isinstance(level, int) and level >= 0, "need positive integer for vibrational level"

        freq = self.frequency
        if freq < MIN_FREQUENCY:
            freq = MIN_FREQUENCY

        # 0.5 * h * c * frequency (c in cm/s bc wavenumbers)
        # 0.5 * (6.626 * 10**-34) * (3 * 10**10) * (6.026 * 10**23) / 4184) = 0.0014305 
        zpe = 0.0014305 * freq
        return zpe * (2 * level + 1)

    def random_displacement(self, energy=None, level=0, method="quasiclassical", max_attempts=1e5):
        """
        Args:
            energy (float): energy of mode (for classical case)
            method (str): "quasiclassical" or "classical"
            level (int): which vibrational level
            max_attempts (int): how many tries you get

        Returns:
            shift
        """
        if method == "quasiclassical":
            min_val = 0
            max_val = self.quantum_distribution_max(level)
            max_x = self.classical_turning_point()

            attempts = 0
            while attempts < max_attempts:
                x = np.random.uniform(-1 * max_x, max_x)
                p = self.quantum_distribution_value(x, level)

                y = np.random.uniform(min_val, max_val)
                if y < p:
                    return x
                else:
                    attempts += 1

            raise ValueError("max_attempts exceeded - can't get a proper initialization for this mode!")
        elif method == "classical":
            assert energy is not None, "need energy for classical displacement"
            min_val = self.classical_distribution_value(0)
            max_x = self.classical_turning_point(energy=energy)
            max_val = self.classical_distribution_value(max_x)

            attempts = 0
            while attempts < max_attempts:
                x = np.random.uniform(-1*max_x, max_x)
                p = self.classical_distribution_value(max_x)

                y = np.random.uniform(min_val, max_val)
                if y < p:
                    return x
                else:
                    attempts += 1
        else:
            raise ValueError(f"invalid method {method} - only ``quasiclassical`` and ``classical`` implemented currently!")

        raise ValueError("Max attempts exceeded!")

    def quantum_distribution_value(self, x, level=0):
        """
        Calculate psi**2 for quantum harmonic oscillator for a given shift in Å.

        Args:
            x (float): shift in Å
            level (int): vibrational level
        """
        assert isinstance(level, int) and level >= 0, "need positive integer for vibrational level"

        freq = self.frequency
        if freq < MIN_FREQUENCY:
            freq = MIN_FREQUENCY

        n = level # brevity is the soul of wit
        H = get_hermite_polynomial(n)

        # following https://github.com/ekwan/Jprogdyn/blob/master/src/main/java/edu/harvard/chemistry/ekwan/Jprogdyn/HarmonicOscillatorDistribution.java, line 109
        # 4 * pi * 3 * 10**8 / (1000 * 10**20 * 6.022 * 10**23 * 6.626 * 10^-34) = 0.000094411, take it or leave it
        omega_term = 9.4411e-5 * self.reduced_mass * freq
        val = math.sqrt(omega_term) * math.exp(-1 * omega_term * math.pi * x ** 2 ) * (H(math.sqrt(omega_term * math.pi) * x) ** 2) / (2 ** n * math.factorial(n))
        return val

    def quantum_distribution_max(self, level=0, num_pts=1e4):
        """
        Returns the maximum value of psi**2 for the quantum harmonic oscillator at a given level.
        """
        assert isinstance(level, int) and level >= 0, "need positive integer for vibrational level"

        if level == 0:
            return self.quantum_distribution_value(0)

        max_x = self.classical_turning_point()

        # there is certainly a better way to do this
        max_p = 0
        for x in np.linspace(0, max_x, int(num_pts)):
            p = self.quantum_distribution_value(x, level)
            if p > max_p:
                max_p = p

        return max_p

    def classical_distribution_value(self, x):
        """
        Returns the value of the classical distribution at the specified ``x`` value.
        """
        max_x = self.classical_turning_point()
        assert (x <= max_x) and (x >= -1*max_x), "x must be in [-max_x, max_x]"
        return 1/(math.pi * math.sqrt(max_x**2 - x**2))

    def classical_turning_point(self, energy=None):
        """
        Returns the maximum allowed shift based on modelling the mode as a classical harmonic oscillator (e.g. the point where potential energy is maximum).

        Args:
            energy (float): energy of mode
            level (int): level to compute energy for quantum harmonic oscillator
        """
        if energy is None:
            energy = self.energy()
        else:
            assert energy > 0, "cannot request turning point for 0 energy!"

        return math.sqrt(2 * energy / self.force_constant)

    def to_string(self):
        ...

    def from_string(self):
        ...
