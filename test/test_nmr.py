import unittest, sys, os, io, copy, math
import numpy as np
import cctk

class TestNMR(unittest.TestCase):
    def load_molecule(self, path="test/static/LSD_custom.out"):
        return cctk.GaussianFile.read_file(path).molecule

    def test_basic(self):
        mol = self.load_molecule()

'''
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
'''

if __name__ == '__main__':
    unittest.main()
