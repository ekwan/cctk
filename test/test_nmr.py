import unittest, sys, os, io, copy, math
import numpy as np
import cctk
from cctk.gaussian_file import GaussianJobType as JobType

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

        tensors = ensemble[:, "shielding_tensors"]
        self.assertEqual(tensors[0][0][0], 198.2259)
        self.assertEqual(tensors[1][0][0], 32.6869)
        for shielding, tensor in zip(shieldings, tensors):
            pred_shielding = np.trace(tensor)/3
            self.assertTrue((shielding-pred_shielding < 0.001))

    def test_nmr2(self):
        # this file contains opt freq followed by Link1 NMR on methane
        gaussian_file = cctk.GaussianFile.read_file("test/static/methane2.out")
        self.assertEqual(len(gaussian_file), 2)
        first_link = gaussian_file[0]
        self.assertLess(abs(first_link.elapsed_time-6.1), 0.001)
        self.assertListEqual(first_link.job_types, [JobType.OPT, JobType.FREQ, JobType.SP])
        ensemble = first_link.ensemble
        energies = list(ensemble[:,"energy"])
        self.assertListEqual(energies, [-40.5169484082, -40.5183831835, -40.5183831835])
        second_link = gaussian_file[1]
        self.assertLess(abs(second_link.elapsed_time-1.9), 0.001)
        ensemble = second_link.ensemble
        shieldings = list(ensemble[-1,"isotropic_shielding"])
        self.assertListEqual(list(shieldings), [192.9242, 31.8851, 31.8851, 31.8851, 31.8851])

        tensors = ensemble[:, "shielding_tensors"]
        for shielding, tensor in zip(shieldings, tensors):
            pred_shielding = np.trace(tensor)/3
            self.assertTrue((shielding-pred_shielding < 0.001))

    def test_nmr3(self):
        # this file contains opt freq / Link1 NMR on ethane then Link1 single point NMR on methane
        gaussian_file = cctk.GaussianFile.read_file("test/static/ethane.out")
        self.assertEqual(len(gaussian_file), 3)

        # opt freq
        first_link = gaussian_file[0]
        self.assertLess(abs(first_link.elapsed_time-11.0), 0.001)
        self.assertListEqual(first_link.job_types, [JobType.OPT, JobType.FREQ, JobType.SP])

        # NMR: ethane
        second_link = gaussian_file[1]
        self.assertLess(abs(second_link.elapsed_time-2.1), 0.001)
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
        self.assertLess(abs(third_link.elapsed_time-1.8), 0.001)
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
        self.assertEqual(gaussian_file.successful_terminations,1)
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

        tensors = ensemble[:, "shielding_tensors"]
        for shielding, tensor in zip(shieldings, tensors):
            pred_shielding = np.trace(tensor)/3
            self.assertTrue((shielding-pred_shielding < 0.001))

    def test_nmr5(self):
        # ensure all acetone files can be read
        for idx in range(1,7):
            if idx==3:
                #### this file is horrible, it should not be read correctly
                with self.assertRaises(ValueError):
                    gaussian_file = cctk.GaussianFile.read_file(f"test/static/acetone-couplings{idx}.out")
            else:
                gaussian_file = cctk.GaussianFile.read_file(f"test/static/acetone-couplings{idx}.out")

    def test_nmr6(self):
        gaussian_file = cctk.GaussianFile.read_file("test/static/acetone-couplings1.out")
        ensemble = gaussian_file.ensemble
        shieldings = ensemble[-1,"isotropic_shielding"]
        expected_shieldings = [165.8515, 30.794, 30.93, 30.9302, -21.4514, -375.1462, 159.4249, 30.6991, 30.6993, 30.8559]
        self.assertTrue((np.abs(shieldings - expected_shieldings) <= 0.0001).all())

        expected_couplings = np.array(\
        [[  0. ,124.1,134.7,134.7, 34.8, -0.8, 15.5, -0.4, -0.4,  4.9],
        [124.1,  0. ,-14.4,-14.4, -3.7, -2.1,  0.6,  0.4,  0.4,  0.7],
        [134.7,-14.4,  0. ,-20.4, -6.9, -1.3,  1.1, -0.6, -1.3, -0.1],
        [134.7,-14.4,-20.4,  0. , -6.9, -1.3,  1.1, -1.2, -0.6, -0.1],
        [ 34.8, -3.7, -6.9, -6.9,  0. , 43.7, 35. , -5.7, -5.6, -6.4],
        [ -0.8, -2.1, -1.3, -1.3, 43.7,  0. , -1.1, -1.9, -1.9, -0.8],
        [ 15.5,  0.6,  1.1,  1.1, 35. , -1.1,  0. ,127.2,127.2,137.5],
        [ -0.4,  0.4, -0.6, -1.2, -5.7, -1.9,127.2,  0. ,-19.1,-14.6],
        [ -0.4,  0.4, -1.3, -0.6, -5.6, -1.9,127.2,-19.1,  0. ,-14.6],
        [  4.9,  0.7, -0.1, -0.1, -6.4, -0.8,137.5,-14.6,-14.6,  0. ]])
        couplings = ensemble[-1,"j_couplings"]
        self.assertTrue(np.any(expected_couplings-couplings < 0.1))

        gaussian_file = cctk.GaussianFile.read_file("test/static/acetone-couplings2.out")
        ensemble = gaussian_file.ensemble
        shieldings = ensemble[-1,"isotropic_shielding"]
        self.assertTrue(shieldings is not None)
        couplings = ensemble[-1,"j_couplings"]
        expected_couplings = np.array(\
        [[  0. ,122.1,132.3,132.3, 39.1,  0.8, 11.6, -0.4, -0.4,  3.5],
        [122.1,  0. ,-14.5,-14.5, -1.8, -0.6,  0.3,  0.4,  0.4,  0.6],
        [132.3,-14.5,  0. ,-19.7, -5.4, -1.5,  0.8, -0.6, -1.1, -0.2],
        [132.3,-14.5,-19.7,  0. , -5.4, -1.5,  0.8, -1.1, -0.6, -0.2],
        [ 39.1, -1.8, -5.4, -5.4,  0. , 42.4, 40.1, -3.6, -3.6, -5.2],
        [  0.8, -0.6, -1.5, -1.5, 42.4,  0. ,  0.3, -1.8, -1.7, -0.1],
        [ 11.6,  0.3,  0.8,  0.8, 40.1,  0.3,  0. ,125.6,125.6,135.2],
        [ -0.4,  0.4, -0.6, -1.1, -3.6, -1.8,125.6,  0. ,-18.6,-14.8],
        [ -0.4,  0.4, -1.1, -0.6, -3.6, -1.7,125.6,-18.6,  0. ,-14.8],
        [  3.5,  0.6, -0.2, -0.2, -5.2, -0.1,135.2,-14.8,-14.8,  0. ]])
        self.assertTrue(np.any(expected_couplings-couplings < 0.1))

        gaussian_file = cctk.GaussianFile.read_file("test/static/acetone-couplings5.out")
        ensemble = gaussian_file[1].ensemble
        shieldings = ensemble[-1,"isotropic_shielding"]
        self.assertTrue(shieldings is not None)
        couplings = ensemble[-1,"j_couplings"]
        expected_couplings = np.array(\
        [[  0. ,126.5,129.7,129.7, 42.3,  0.2, 13.5, -0.4, -0.4,  4.5],
        [126.5,  0. ,-14.4,-14.4, -0. , -0.9,  0.2,  0.5,  0.5,  0.8],
        [129.7,-14.4,  0. ,-21.2, -4.7, -1.4,  0.9, -0. , -1. , -0. ],
        [129.7,-14.4,-21.2,  0. , -4.7, -1.4,  0.9, -1. , -0. , -0. ],
        [ 42.3, -0. , -4.7, -4.7,  0. , 37.4, 44.1, -2.6, -2.6, -4.8],
        [  0.2, -0.9, -1.4, -1.4, 37.4,  0. , -0.2, -1.7, -1.7, -0.3],
        [ 13.5,  0.2,  0.9,  0.9, 44.1, -0.2,  0. ,125.6,125.6,135.6],
        [ -0.4,  0.5, -0. , -1. , -2.6, -1.7,125.6,  0. ,-20.6,-14.3],
        [ -0.4,  0.5, -1. , -0. , -2.6, -1.7,125.6,-20.6,  0. ,-14.3],
        [  4.5,  0.8, -0. , -0. , -4.8, -0.3,135.6,-14.3,-14.3,  0. ]])
        self.assertTrue(np.any(expected_couplings-couplings < 0.1))

        gaussian_file = cctk.GaussianFile.read_file("test/static/acetone-couplings6.out")
        ensemble = gaussian_file[1].ensemble
        shieldings = ensemble[-1,"isotropic_shielding"]
        self.assertTrue(shieldings is not None)
        couplings = ensemble[-1,"j_couplings"]
        expected_couplings = np.array(\
        [[  0. ,130.6,130.7,130.8, 38.4, -1.6, 19.1, -0.4, -0.4,  6. ],
        [130.6,  0. ,-14.2,-14.2, -1.6, -2.5,  0.5,  0.4,  0.4,  1. ],
        [130.7,-14.2,  0. ,-22. , -6.3, -1.2,  1.2,  0. , -1. ,  0.1],
        [130.8,-14.2,-22. ,  0. , -6.3, -1.2,  1.2, -1.1,  0. ,  0.1],
        [ 38.4, -1.6, -6.3, -6.3,  0. , 33.8, 39.1, -4.8, -4.8, -5.9],
        [ -1.6, -2.5, -1.2, -1.2, 33.8,  0. , -1.8, -1.9, -1.9, -1. ],
        [ 19.1,  0.5,  1.2,  1.2, 39.1, -1.8,  0. ,126.6,126.6,138.3],
        [ -0.4,  0.4,  0. , -1.1, -4.8, -1.9,126.6,  0. ,-21.3,-14. ],
        [ -0.4,  0.4, -1. ,  0. , -4.8, -1.9,126.6,-21.3,  0. ,-14. ],
        [  6. ,  1. ,  0.1,  0.1, -5.9, -1. ,138.3,-14. ,-14. ,  0. ]])
        self.assertTrue(np.any(expected_couplings-couplings < 0.1))

        #with np.printoptions(precision=1, suppress=True):
        #    print(np.array2string(couplings, separator=","))

    def test_nmr_solvated(self):
        f1 = cctk.GaussianFile.read_file("test/static/ibuprofen_solvated.out")
        f2 = cctk.GaussianFile.read_file("test/static/ibuprofen_solvated2.out")

if __name__ == '__main__':
    unittest.main()
