import unittest, sys, os, io, copy
import numpy as np
import cctk

# python -m unittest test.test_align.TestAlign
class TestAlign(unittest.TestCase):
    def test_align(self):
        path = "test/static/gaussian_file.out"
        gaussian_file = cctk.GaussianFile.read_file(path)
        conformational_ensemble = gaussian_file.molecules
        m1 = conformational_ensemble[0]
        m2 = conformational_ensemble[-1]
        RMSD = cctk.helper_functions.compute_RMSD(m1,m2)
        delta = abs(0.0006419131435567976 - RMSD)
        self.assertLess(delta, 0.0001)

if __name__ == '__main__':
    unittest.main()
