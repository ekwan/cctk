import unittest, sys, os, io, copy
import numpy as np
import cctk
import glob as glob

# python -m unittest test.test_ensemble2.TestEnsemble2
class TestEnsemble2(unittest.TestCase):
    def test_ensemble(self):
        path = "test/static/phenylpropane*.out"
        conformational_ensemble = cctk.ConformationalEnsemble()
        for filename in sorted(glob.glob(path)):
            gaussian_file = cctk.GaussianFile.read_file(filename)
            e = gaussian_file.molecules
            m = e.molecules[-1]
            p = e[m]
            conformational_ensemble.add_molecule(m,p)
        print(conformational_ensemble[0])
        print(conformational_ensemble[0:2])

if __name__ == '__main__':
    unittest.main()
