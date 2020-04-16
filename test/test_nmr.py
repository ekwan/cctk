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
        ensemble = gaussian_file.ensemble
        molecule = ensemble.molecules[-1]
        properties_dict = ensemble.get_properties_dict(molecule)
        energy = properties_dict["energy"]
        self.assertEqual(energy, -40.5169484082)
        self.assertEqual(ensemble[-1,"energy"], -40.5169484082)
        shieldings = properties_dict["isotropic_shielding"]
        self.assertListEqual(list(shieldings), [198.2259, 32.6869, 32.6869, 32.6869, 32.6869])
        shieldings = ensemble[:,"isotropic_shielding"]
        self.assertListEqual(list(shieldings), [198.2259, 32.6869, 32.6869, 32.6869, 32.6869])

    def test_nmr2(self):
        # this file contains opt freq followed by Link1 NMR on methane
        gaussian_file = cctk.GaussianFile.read_file("test/static/methane2.out")
        self.assertEqual(len(gaussian_file), 2)
        first_link = gaussian_file[0]
        self.assertListEqual(first_link.job_types, [JobType.OPT, JobType.FREQ, JobType.SP])
        ensemble = first_link.ensemble
        energies = list(ensemble[:,"energy"])
        self.assertListEqual(energies, [-40.5169484082, -40.5183831835, -40.5183831835])
        second_link = gaussian_file[1]
        ensemble = second_link.ensemble
        shieldings = list(ensemble[-1,"isotropic_shielding"])
        self.assertListEqual(list(shieldings), [192.9242, 31.8851, 31.8851, 31.8851, 31.8851])

    def test_nmr3(self):
        # this file contains opt freq / Link1 NMR on ethane then Link1 single point NMR on methane
        gaussian_file = cctk.GaussianFile.read_file("test/static/ethane.out")
        self.assertEqual(len(gaussian_file), 3)

        # opt freq
        first_link = gaussian_file[0]
        self.assertListEqual(first_link.job_types, [JobType.OPT, JobType.FREQ, JobType.SP])

        # NMR: ethane
        second_link = gaussian_file[1]
        ensemble = second_link.ensemble
        molecule = ensemble[-1]
        #shifts = list(ensemble[molecule]["isotropic_shielding"])
        shifts = ensemble[:,"isotropic_shielding"]  # list
        self.assertListEqual(list(shifts), [180.3673, 31.2068, 31.207, 31.2068, 180.3673, 31.2068, 31.207, 31.2068])
        scaled_shifts, shift_labels = cctk.helper_functions.scale_nmr_shifts(ensemble,
                                      symmetrical_atom_numbers=[[1,5],[2,3,4,6,7,8]], scaling_factors="default")
        self.assertTrue((np.abs(scaled_shifts[0] - np.array([0.42845589, 0.06087379])) <= 0.00001).all())

        # NMR: methane
        third_link = gaussian_file[2]
        self.assertListEqual(third_link.job_types, [JobType.NMR, JobType.SP])
        ensemble = third_link.ensemble
        #molecule = ensemble[-1]
        #shifts = list(ensemble[molecule]["isotropic_shielding"])
        shieldings = ensemble[:,"isotropic_shielding"]
        self.assertListEqual(list(shieldings), [198.2259, 32.6869, 32.6869, 32.6869, 32.6869])
        scaled_shifts, shift_labels = cctk.helper_functions.scale_nmr_shifts(ensemble,
                                      symmetrical_atom_numbers=None, scaling_factors="default")

    def test_nmr4(self):
        # tests code for scaling NMR shieldings on a more complicated molecule
        gaussian_file = cctk.GaussianFile.read_file("test/static/LSD_custom.out")
        ensemble = gaussian_file.ensemble
        shieldings = ensemble[:,"isotropic_shielding"]
        scaled_shifts, shift_labels = cctk.helper_functions.scale_nmr_shifts(ensemble,
                                      symmetrical_atom_numbers=[[37,38,39],[32,33,34]], scaling_factors="default")
        expected_shifts = [6.52352,6.6285,6.51045,6.53005,6.22303,2.11021,2.7025,2.73022,2.38541,2.35172,3.1467,5.82979,
                           3.29202,1.92326,114.924,98.3836,107.641,94.3333,104.421,109.795,95.1041,112.168,121.346,
                           45.4898,14.1014,26.7028,36.3779,29.4323,104.708,155.804,38.0661,109.579,22.7099]
        expected_shifts = np.asarray(expected_shifts)
        self.assertTrue((np.abs(scaled_shifts[0] - expected_shifts) <= 0.001).all())
        #print(shift_labels)

if __name__ == '__main__':
    unittest.main()
