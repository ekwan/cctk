import unittest, sys, os, io, copy, math, shutil
import numpy as np

import cctk
import cctk.optimize as opt

class TestMolecule(unittest.TestCase):
    def load_molecule(self, path="test/static/test_peptide.xyz"):
        return cctk.XYZFile.read_file(path).get_molecule()

    def test_basic(self):
        mol = self.load_molecule()
#        e1 = mol.compute_energy()
#        print(e1)

        mol2, e2 = opt.optimize_molecule(mol, nprocs=4, return_energy=True)
        self.assertTrue(isinstance(mol2, cctk.Molecule))
        self.assertTrue(e2 + 66.394 < 0.1)

        mol.optimize(nprocs=4) #in-place

        for x1, x2 in zip(np.ravel(mol2.geometry), np.ravel(mol.geometry)):
            self.assertTrue(abs(float(x1)-float(x2)) < 0.1)

        self.assertTrue(mol.compute_energy() + 66.474523114714 < 0.00001)

    def skip_test_csearch(self):
        mol = cctk.GaussianFile.read_file("test/static/L-Ala.gjf").get_molecule()

        if opt.installed("crest") is not None:
            conformers = mol.csearch()
            self.assertEqual(len(conformers), 36)
            self.assertTrue(isinstance(conformers, cctk.ConformationalEnsemble))
        else:
            pass


if __name__ == '__main__':
    unittest.main()
