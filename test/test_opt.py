import unittest, sys, os, io, copy, math
import numpy as np

import cctk
import cctk.optimize as opt

class TestMolecule(unittest.TestCase):
    def load_molecule(self, path="test/static/test_peptide.xyz"):
        return cctk.XYZFile.read_file(path).molecule

    def test_basic(self):
        mol = self.load_molecule()
        mol2 = opt.optimize_molecule(mol, nprocs=4)
        self.assertTrue(isinstance(mol2, cctk.Molecule))

        mol.optimize(nprocs=4) #in-place

        for x1, x2 in zip(np.ravel(mol2.geometry), np.ravel(mol.geometry)):
            self.assertTrue(abs(float(x1)-float(x2)) < 0.1)

if __name__ == '__main__':
    unittest.main()
