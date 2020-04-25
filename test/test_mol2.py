import unittest, sys, os, io, copy
import numpy as np
import cctk
import time

class TestMOL2(unittest.TestCase):
    def test_read(self):
        path = "test/static/dodecane.mol2"
        file = cctk.MOL2File.read_file(path)
        self.assertTrue(isinstance(file, cctk.MOL2File))

        ensemble = file.ensemble
        self.assertTrue(isinstance(ensemble, cctk.ConformationalEnsemble))
        self.assertEqual(len(ensemble), 1)

        mol = ensemble.molecules[0]
        self.assertTrue(isinstance(mol, cctk.Molecule))
        self.assertEqual(len(mol.atomic_numbers), 38)
        self.assertEqual(len(mol.geometry), 38)
        self.assertEqual(mol.get_bond_order(1,2), 1)

    def test_bulk_read(self):
        path = "test/static/dodecane-csearch.mol2"
        file = cctk.MOL2File.read_file(path)
        self.assertTrue(isinstance(file, cctk.MOL2File))

        ensemble = file.ensemble
        self.assertTrue(isinstance(ensemble, cctk.ConformationalEnsemble))
        self.assertEqual(len(ensemble), 597)
        for mol in ensemble.molecules:
            self.assertEqual(len(mol.atomic_numbers), 38)
            self.assertEqual(len(mol.geometry), 38)
            self.assertEqual(mol.get_bond_order(1,2), 1)

    def test_bulk_read_big(self):
        with open("test/static/dodecane-csearch.mol2", "r") as source_file:
            lines = source_file.readlines()
            with open("test/static/dodecane-csearch-big.mol2", "w") as destination_file:
                for i in range(20):
                    destination_file.writelines(lines)
        path = "test/static/dodecane-csearch-big.mol2"
        start = time.time()
        file = cctk.MOL2File.read_file(path, print_status_messages=False, contains_conformers=True)
        end = time.time()
        delta = end-start
        #print(f"elapsed time: {delta:.3f} s")
        #print(file.ensemble)
        os.remove(path)

    def test_regular_ensemble(self):
        path = "test/static/dodecane2.mol2"
        file = cctk.MOL2File.read_file(path)
        self.assertTrue(isinstance(file, cctk.MOL2File))
        ensemble = file.ensemble
        self.assertTrue(isinstance(ensemble, cctk.Ensemble))
        self.assertEqual(len(ensemble), 2)
        self.assertFalse(isinstance(ensemble, cctk.ConformationalEnsemble))

    def test_write(self):
        path = "test/static/adamantane.mol2"
        file = cctk.MOL2File.read_file(path)
        new_path = "test/static/new_Ad.mol2"

        file.write_file(new_path, title="title")

        with open(path) as old:
            with open(new_path) as new:
                self.assertListEqual(
                    list(new),
                    list(old)
                )

        os.remove(new_path)

if __name__ == '__main__':
    unittest.main()
