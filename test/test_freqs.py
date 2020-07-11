import cctk
import unittest
import numpy as np

import cctk.quasiclassical as qc

if __name__ == '__main__':
    unittest.main()

class TestFrequencies(unittest.TestCase):
    def test_read(self):
        path = "test/static/gaussian_file.out"
        file = cctk.GaussianFile.read_file(path)

        mol = file.get_molecule()
        self.assertEqual(len(mol.vibrational_modes), 87)
        self.assertEqual(mol.vibrational_modes[-1].frequency, 3614.2603)

        mol2 = qc.get_quasiclassical_vibration(mol)

    def test_qho(self):
        from asciichartpy import plot
        path = "test/static/gaussian_file.out"
        file = cctk.GaussianFile.read_file(path)
        mol = file.get_molecule()

        mode = mol.vibrational_modes[-1]

        for level in range(10):
            tp = mode.classical_turning_point(level=level)
            print(f"LEVEL {level} ({mode.energy(level):.2f} kcal.mol, turning point = {tp:.2f} Å)")
            tp = 2
            print(plot([mode.quantum_distribution_value(x, level=level) for x in np.linspace(-1 * tp, tp, 100)], {"height": 8}))

