import unittest, sys, os, io, copy, math
import numpy as np
import cctk
from cctk.gaussian_file import JobType

# run from cctk root with
# python -m unittest test.test_nmr.TestNMR
class TestNMR(unittest.TestCase):

    def test_nmr1(self):
        # this file has a single point NMR calculation on methane
        gaussian_file = cctk.GaussianFile.read_file("test/static/methane.out")
        ensemble = gaussian_file.molecules
        molecule = ensemble[-1]
        self.assertListEqual(list(molecule.nmr_isotropic), [198.2259, 32.6869, 32.6869, 32.6869, 32.6869])

    def test_nmr2(self):
        # this file contains opt freq followed by Link1 NMR on methane
        gaussian_file = cctk.GaussianFile.read_file("test/static/methane2.out")
        self.assertEqual(len(gaussian_file), 2)
        first_link = gaussian_file[0]
        self.assertListEqual(first_link.job_types, [JobType.OPT, JobType.FREQ])
        second_link = gaussian_file[1]
        ensemble = second_link.molecules
        molecule = ensemble[-1]
        self.assertListEqual(list(molecule.nmr_isotropic), [192.9242, 31.8851, 31.8851, 31.8851, 31.8851])

    def test_nmr3(self):
        # this file contains opt freq / Link1 NMR on ethane then Link1 single point NMR on methane
        gaussian_file = cctk.GaussianFile.read_file("test/static/ethane.out")
        self.assertEqual(len(gaussian_file), 3)
        first_link = gaussian_file[0]
        self.assertListEqual(first_link.job_types, [JobType.OPT, JobType.FREQ])
        second_link = gaussian_file[1]
        self.assertListEqual(second_link.job_types, [JobType.NMR])
        ensemble = second_link.molecules
        molecule = ensemble[-1]
        self.assertListEqual(list(molecule.nmr_isotropic), [180.3673, 31.2068, 31.207, 31.2068, 180.3673, 31.2068, 31.207, 31.2068])
        third_link = gaussian_file[2]
        self.assertListEqual(third_link.job_types, [JobType.NMR])
        ensemble = third_link.molecules
        molecule = ensemble[-1]
        self.assertListEqual(list(molecule.nmr_isotropic), [198.2259, 32.6869, 32.6869, 32.6869, 32.6869])

'''
    def load_molecule(self, path="test/static/LSD_custom.out"):
        return cctk.GaussianFile.read_file(path).get_molecule()

    def test_basic(self):
        mol = self.load_molecule()

        self.assertTrue(isinstance(mol.nmr_isotropic, cctk.OneIndexedArray))
        self.assertEqual(len(mol.nmr_isotropic), len(mol.atomic_numbers))

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
