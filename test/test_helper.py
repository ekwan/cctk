import unittest, sys, os, io, copy
import numpy as np
import cctk

import cctk.helper_functions as helper

class TestHelper(unittest.TestCase):
    def test_get_number(self):
        self.assertEqual(helper.get_number("CL"), 17)
        self.assertEqual(helper.get_number("cl"), 17)
        self.assertEqual(helper.get_number("Cl"), 17)

        self.assertEqual(helper.get_number("Bq"), 0)
        self.assertEqual(helper.get_number("U"), 92)

    def test_get_symbol(self):
        self.assertEqual(helper.get_symbol(0), "Bq")
        self.assertEqual(helper.get_symbol(62), "Sm")
        self.assertEqual(helper.get_symbol(46), "Pd")
        self.assertEqual(helper.get_symbol(6), "C")

    def test_isotope(self):
        m, w = helper.get_isotopic_distribution(1)
        self.assertListEqual(list(m), [1.007825, 2.014102, 3.016049])
        self.assertListEqual(list(w), [0.999885, 0.000115, 0.0])

    def test_free_energy(self):
        file = cctk.GaussianFile.read_file("test/static/diiron_complex.out")
        ensemble = file.ensemble
        molecule = ensemble.molecules[-1]
        properties_dict = ensemble.get_properties_dict(molecule)
        #print(properties_dict)
        free_energy = properties_dict["gibbs_free_energy"]
        corrected_free_energy = properties_dict["quasiharmonic_gibbs_free_energy"]
        frequencies = properties_dict["frequencies"]
        self.assertLess(abs(free_energy+1975.998622), 0.00001)
        self.assertLess(abs(corrected_free_energy+1975.9949313424), 0.00001)
        corrected_free_energy = cctk.helper_functions.get_corrected_free_energy(free_energy, frequencies,
                                                                                frequency_cutoff=50.0, temperature=298.15)
        self.assertLess(abs(corrected_free_energy+1975.997664138452), 0.00001)
        #delta = (-free_energy+corrected_free_energy)*627.509469
        #print(delta)

if __name__ == '__main__':
    unittest.main()
