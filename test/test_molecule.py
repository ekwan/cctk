import unittest, sys, os, io, copy, math
import numpy as np
import cctk

class TestMolecule(unittest.TestCase):
    def load_molecule(self, path="static/test_peptide.xyz"):
        return cctk.XYZFile.read_file(path).molecule

    def test_basic(self):
        mol = self.load_molecule()
        self.assertListEqual(mol.atomic_numbers.tolist(), [7, 1, 6, 1, 6, 6, 1, 8, 7, 1, 6, 1, 6, 6, 1, 8, 8, 6, 1, 1, 1, 9, 9, 9, 9, 6, 8, 6, 1, 1, 1])

        mol.assign_connectivity()
        edges = [(1, 2), (1, 3), (1, 26), (3, 4), (3, 5), (3, 6), (5, 7), (5, 24), (5, 25), (6, 8), (6, 9), (9, 10), (9, 11), (11, 12), (11, 13), (11, 14), (13, 15), (13, 22), (13, 23), (14, 16), (14, 17), (17, 18), (18, 19), (18, 20), (18, 21), (26, 27), (26, 28), (28, 29), (28, 30), (28, 31)]
        self.assertListEqual(list(mol.bonds.edges()), edges)

    def test_distance(self):
        mol = self.load_molecule()
        self.assertTrue(isinstance(mol, cctk.Molecule))

        self.assertEqual(int(round(mol.get_distance(1,2)*10)), 10)
        self.assertEqual(int(round(mol.get_distance(1,3)*10)), 14)
        self.assertEqual(int(round(mol.get_distance(1,9)*10)), 38)

        mol.set_distance(1, 2, 2.00)

        self.assertEqual(int(round(mol.get_distance(1,2)*10)), 20)
        self.assertEqual(int(round(mol.get_distance(1,3)*10)), 14)
        self.assertEqual(int(round(mol.get_distance(1,9)*10)), 38)

        self.assertTrue(mol.check_for_conflicts())
        mol.set_distance(1, 2, 0.01)
        self.assertRaises(ValueError, mol.check_for_conflicts)

    def test_angle(self):
        mol = self.load_molecule()
        self.assertTrue(isinstance(mol, cctk.Molecule))

        self.assertEqual(int(round(mol.get_angle(1,3,5))), 111)
        self.assertEqual(int(round(mol.get_angle(3,5,7))), 110)
        self.assertEqual(int(round(mol.get_angle(5,7,9))), 64)

        mol.set_angle(1, 3, 5, 120)

        self.assertEqual(int(round(mol.get_angle(1,3,5))), 120)

    def test_dihedral(self):
        mol = self.load_molecule()
        self.assertTrue(isinstance(mol, cctk.Molecule))

        self.assertEqual(int(round(mol.get_dihedral(1,3,5,7))), 60)
        self.assertEqual(int(round(mol.get_dihedral(16,14,17,18))), 11)
        self.assertEqual(int(round(mol.get_dihedral(31,28,1,2))), 148)

        theta = [1, 20, 89, 66, 180, 215, 333]
        for t in theta:
            mol.set_dihedral(1, 3, 5, 7, t)
            self.assertEqual(int(round(mol.get_dihedral(1,3,5,7))), t)

    def test_fragment(self):
        mol = self.load_molecule()
        mol.assign_connectivity()
        (frag1, frag2) = mol._get_bond_fragments(3, 5)
        self.assertEqual(len(frag1), 27)
        self.assertEqual(len(frag2), 4)

        self.assertEqual(len(mol._get_fragment_containing(5)), 31)
        mol.remove_bond(3,5)
        self.assertEqual(len(mol._get_fragment_containing(5)), 4)
        self.assertFalse(mol.are_connected(3,5))
        mol.add_bond(3,5)
        self.assertEqual(len(mol._get_fragment_containing(5)), 31)
        self.assertTrue(mol.are_connected(3,5))

    def test_add_atoms(self):
        mol = cctk.Molecule(np.array([2], dtype=np.int8), [[0, 0, 0]])
        self.assertEqual(mol.num_atoms(), 1)

        mol.add_atom("He", [1, 0, 0])
        self.assertListEqual(mol.atomic_numbers.tolist(), [2, 2])
        self.assertEqual(mol.num_atoms(), 2)

        mol.add_atom("Ar", [3, 0, 0])
        self.assertEqual(mol.num_atoms(), 3)

        mol.add_atom_at_centroid("He", [2, 3])
        self.assertEqual(mol.num_atoms(), 4)
        self.assertListEqual(list(mol.get_vector(4)), [2, 0, 0])

    def test_mass_spec(self):
        mol = cctk.Molecule(np.array([12], dtype=np.int8), [[0, 0, 0]])
        masses, weights = mol.calculate_mass_spectrum()
        self.assertListEqual(list(masses), [23.])
        self.assertListEqual(list(weights), [1.])

    def test_selection(self):
        mol = self.load_molecule()
        self.assertListEqual(mol.get_heavy_atoms(), [1, 3, 5, 6, 8, 9, 11, 13, 14, 16, 17, 18, 22, 23, 24, 25, 26, 27, 28])
        self.assertListEqual(mol.get_atoms_by_symbol("F"), [22, 23, 24, 25])

    def test_translate(self):
        mol = cctk.Molecule(np.array([12], dtype=np.int8), [[0, 0, 0]])

        v = np.array([1.5234,1.231234,-1.77777])
        mol = mol.translate_molecule(v)

        self.assertListEqual(mol.geometry.tolist()[0], list(v))
        self.assertTrue(isinstance(mol.geometry, cctk.OneIndexedArray))

        mol2 = self.load_molecule()
        v2 = np.zeros(shape=3)

        mol2_shift = mol2.translate_molecule(v2)
        self.assertListEqual(mol2.geometry.tolist()[0], mol2_shift.geometry.tolist()[0])

if __name__ == '__main__':
    unittest.main()
