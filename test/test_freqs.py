import cctk
import unittest
import numpy as np

import cctk.quasiclassical as qc

if __name__ == '__main__':
    unittest.main()

class TestFrequencies(unittest.TestCase):
    def test_read(self):
        path = "test/static/methane_normal.out"
        file = cctk.GaussianFile.read_file(path)

        mol = file.get_molecule()
        self.assertEqual(len(mol.vibrational_modes), 9)
        self.assertEqual(mol.vibrational_modes[-1].frequency, 3119.6807)
        self.assertEqual(mol.vibrational_modes[-1].intensity, 26.1449)

        path = "test/static/methane_hpmodes.out"
        file = cctk.GaussianFile.read_file(path)
        mol = file.get_molecule()
        self.assertEqual(len(mol.vibrational_modes), 9)
        self.assertEqual(mol.vibrational_modes[-1].frequency, 3121.6678)
        self.assertEqual(mol.vibrational_modes[-1].intensity, 26.1383)

    def draw_qho(self):
        from asciichartpy import plot

        path = "test/static/methane_normal.out"
        file = cctk.GaussianFile.read_file(path)
        mol = file.get_molecule()

        mode = mol.vibrational_modes[-1]

        for level in range(10):
            tp = mode.classical_turning_point(level=level)
            print(f"LEVEL {level} ({mode.energy(level):.2f} kcal.mol, turning point = {tp:.2f} Å)")
            tp = 5
            print(plot([mode.quantum_distribution_value(x, level=level) for x in np.linspace(-1 * tp, tp, 100)], {"height": 10}))

    def test_perturb(self):
        path = "test/static/methane_hpmodes.out"
        file = cctk.GaussianFile.read_file(path)

        mol = file.get_molecule()

        energies = list()
        for i in range(1000):
            _, e, _, _ = qc.get_quasiclassical_perturbation(mol)
            energies.append(e)

        energies = np.array(energies)
#        from asciichartpy import plot
#        hist, edges = np.histogram(energies, bins=120)
#        print(plot(hist, {"height": 10}))

        self.assertTrue(abs(np.mean(energies) - 9.4) < 0.5)
        self.assertTrue(abs(np.std(energies) - 3) < 0.5)

    def test_perturb_water(self):
        path = "test/static/h2o_hpmodes.out"
        file = cctk.GaussianFile.read_file(path)

        mol = file.get_molecule()
        mol2, e, _, text = qc.get_quasiclassical_perturbation(mol)
        self.assertTrue(isinstance(mol2, cctk.Molecule))

        mo = {
            1: {"velocity": "zero"},
            2: {"velocity": "zero"},
            3: {"velocity": "zero"},
        }

        mol3, e, _, _, v = qc.get_quasiclassical_perturbation(mol, return_velocities=True, mode_options=mo)
        self.assertTrue(isinstance(mol3, cctk.Molecule))
        self.assertFalse(np.any(v)) # all should be zero, AKA False

        mo = {
            1: {"velocity": "positive", "displacement": False},
            2: {"velocity": "positive", "displacement": False},
            3: {"velocity": "positive", "displacement": False},
        }
        mol4, e, te, text, v = qc.get_quasiclassical_perturbation(mol, return_velocities=True, mode_options=mo)
        self.assertTrue(te - 13.28839457 < 0.00001)

        mol5, e, te, text, v = qc.get_quasiclassical_perturbation(mol, return_velocities=True, which="classical")

    def test_final_structure(self):
        path1 = "test/static/methane_perturbed.gjf"
        path2 = "test/static/methane_perturbed_key.gjf"

        e = cctk.ConformationalEnsemble()

        e.add_molecule(cctk.GaussianFile.read_file(path1).get_molecule().assign_connectivity())
        e.add_molecule(cctk.GaussianFile.read_file(path2).get_molecule().assign_connectivity().renumber_to_match(e.molecules[0]))

        e2, before, after = e.align(comparison_atoms="all", compute_RMSD=True)

        self.assertTrue(after[1] < 0.005)

    def test_imaginary(self):
        path = "test/static/eliminationTS.out"
        file = cctk.GaussianFile.read_file(path)
        mol = file.get_molecule()

        self.assertListEqual(mol.atoms_moving_in_imaginary(), [48,2,51])

    def test_boltzmann(self):
        energies = np.zeros(shape=10000)
        for i in range(10000):
            energies[i] = cctk.quasiclassical.random_boltzmann_energy(298)
        mean_energy = np.average(energies)

        # should be 0.2965 = 0.5 kT
        self.assertTrue(abs(0.2965 - mean_energy) < 0.2)

