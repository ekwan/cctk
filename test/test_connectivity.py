import unittest, sys, os, io, copy
import numpy as np
import cctk

import cctk.helper_functions as helper

class TestConnectivity(unittest.TestCase):
    def test_autoassign(self):
        mol = cctk.GaussianFile.read_file("test/static/L-Ala.gjf").get_molecule()
        mol.assign_connectivity(0.1)
        self.assertEqual(len(mol.bonds.edges()), 12)

        mol = cctk.GaussianFile.read_file("test/static/renumber_0.gjf").get_molecule()
        mol.assign_connectivity(0.1)
        self.assertEqual(len(mol.bonds.edges()), 31)
if __name__ == '__main__':
    unittest.main()
