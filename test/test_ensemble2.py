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
            p = e[m].properties_list()[0]
            conformational_ensemble.add_molecule(m,p)
        print(conformational_ensemble.molecules[0])
        print(conformational_ensemble.molecules[0:2])
        print(conformational_ensemble.molecules[[0,1,2]])

    def test_ensemble_indexing(self):
        path = "test/static/gaussian_file.out"
        file = cctk.GaussianFile.read_file(path)
        mols = file.molecules
        self.assertTrue(isinstance(mols, cctk.ConformationalEnsemble))

        self.assertEqual(len(mols), 3)
        self.assertTrue(isinstance(mols[0], cctk.ConformationalEnsemble))

        self.assertListEqual(mols.get_property(None, "energy"), [-1159.56782625, -1159.56782622, -1159.56782622])
        self.assertListEqual(mols.get_property(None, "enthalpy"), [None, None, -1159.314817])

        self.assertListEqual(mols[:,"energy"], [-1159.56782625, -1159.56782622, -1159.56782622])
        self.assertListEqual(mols[:,"enthalpy"], [None, None, -1159.314817])

        self.assertEqual(mols[-1,"energy"], -1159.56782622)
        self.assertEqual(mols[-1,"enthalpy"], -1159.314817)
        self.assertEqual(mols[2, "energy"], -1159.56782622)
        self.assertEqual(mols[2, "enthalpy"], -1159.314817)

if __name__ == '__main__':
    unittest.main()
