import unittest
import cctk

class TestMAE(unittest.TestCase):
    def test_read(self):
        path = "test/static/dodecane_csearch-out.mae"
        (file, pnames, pvals) = cctk.MAEFile.read_file(path)
        self.assertTrue(isinstance(file, cctk.MAEFile))
        self.assertEqual(len(pnames), 597)
        self.assertEqual(len(pvals), 597)

        ensemble = file.ensemble
        self.assertTrue(isinstance(ensemble, cctk.ConformationalEnsemble))
        self.assertEqual(len(ensemble), 597)
        for mol in ensemble.molecules:
            self.assertEqual(len(mol.atomic_numbers), 38)
            self.assertEqual(len(mol.geometry), 38)
            self.assertEqual(mol.get_bond_order(1,2), 1)

if __name__ == '__main__':
    unittest.main()
