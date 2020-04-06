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
        properties = ensemble[molecule]
        energy = properties["energy"]
        self.assertEqual(energy, -40.5169484082)
        shieldings = properties["isotropic_shielding"]
        self.assertListEqual(list(shieldings), [198.2259, 32.6869, 32.6869, 32.6869, 32.6869])
        shieldings = ensemble[:,"isotropic_shielding"]
        self.assertListEqual(list(shieldings[0]), [198.2259, 32.6869, 32.6869, 32.6869, 32.6869])

    def test_nmr2(self):
        # this file contains opt freq followed by Link1 NMR on methane
        gaussian_file = cctk.GaussianFile.read_file("test/static/methane2.out")
        self.assertEqual(len(gaussian_file), 2)
        first_link = gaussian_file[0]
        self.assertListEqual(first_link.job_types, [JobType.OPT, JobType.FREQ, JobType.SP])
        ensemble = first_link.molecules
        #energies = [ ensemble[molecule]["energy"] for molecule in ensemble.molecules() ]
        energies = list(ensemble[:,"energy"])
        self.assertListEqual(energies, [-40.5169484082, -40.5183831835, -40.5183831835])
        #for molecule,properties in ensemble:
        #    print(molecule)
        #    print(properties)
        second_link = gaussian_file[1]
        ensemble = second_link.molecules
        #last_molecule = ensemble[-1]
        shifts = list(ensemble[-1,"isotropic_shielding"])
        self.assertListEqual(shifts, [192.9242, 31.8851, 31.8851, 31.8851, 31.8851])

    def test_nmr3(self):
        # this file contains opt freq / Link1 NMR on ethane then Link1 single point NMR on methane
        gaussian_file = cctk.GaussianFile.read_file("test/static/ethane.out")
        self.assertEqual(len(gaussian_file), 3)
        first_link = gaussian_file[0]
        self.assertListEqual(first_link.job_types, [JobType.OPT, JobType.FREQ, JobType.SP])
        second_link = gaussian_file[1]
        self.assertListEqual(second_link.job_types, [JobType.NMR, JobType.SP])
        ensemble = second_link.molecules
        molecule = ensemble[-1]
        #shifts = list(ensemble[molecule]["isotropic_shielding"])
        shifts = ensemble[:,"isotropic_shielding"]
        self.assertListEqual(list(shifts[0]), [180.3673, 31.2068, 31.207, 31.2068, 180.3673, 31.2068, 31.207, 31.2068])
        third_link = gaussian_file[2]
        self.assertListEqual(third_link.job_types, [JobType.NMR, JobType.SP])
        ensemble = third_link.molecules
        #molecule = ensemble[-1]
        #shifts = list(ensemble[molecule]["isotropic_shielding"])
        shifts = ensemble[:,"isotropic_shielding"]
        self.assertListEqual(list(shifts[0]), [198.2259, 32.6869, 32.6869, 32.6869, 32.6869])

if __name__ == '__main__':
    unittest.main()
