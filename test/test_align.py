import unittest, sys, os, io, copy
import numpy as np
import cctk

# python -m unittest test.test_align.TestAlign
class TestAlign(unittest.TestCase):
    def test_RMSD(self):
        path = "test/static/gaussian_file.out"
        gaussian_file = cctk.GaussianFile.read_file(path)
        conformational_ensemble = gaussian_file.molecules
        m1 = conformational_ensemble[0]
        m2 = conformational_ensemble[-1]
        RMSD = cctk.helper_functions.compute_RMSD(m1,m2)
        delta = abs(0.0006419131435567976 - RMSD)
        self.assertLess(delta, 0.0001)

    def test_align(self):
        path = "test/static/gaussian_file.out"
        gaussian_file = cctk.GaussianFile.read_file(path)
        conformational_ensemble = gaussian_file.molecules
        aligned_ensemble, before_RMSD, after_RMSD = conformational_ensemble.align(to_geometry=0, comparison_atoms="heavy", compute_RMSD=True)
        for i,j in zip(before_RMSD, after_RMSD):
            print(f"{i:6.2f}   {j:6.2f}")

if __name__ == '__main__':
    unittest.main()
