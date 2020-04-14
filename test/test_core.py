import unittest, sys, os, io, copy
import numpy as np
import cctk

class TestXYZ(unittest.TestCase):
    def test_readfile(self):
        path = "test/static/test_peptide.xyz"
        file = cctk.XYZFile.read_file(path)
        self.assertEqual(file.title, "peptide example")

        mol = file.molecule
        self.assertTrue(isinstance(mol, cctk.Molecule))
        self.assertEqual(mol.num_atoms(), 31)
        self.assertTrue(mol.check_for_conflicts())

    def test_writefile(self):
        path = "test/static/test_peptide.xyz"
        new_path = "test/static/test_peptide_copy.xyz"

        file = cctk.XYZFile.read_file(path)
        file.write_file(new_path)

        with open(path) as old:
            with open(new_path) as new:
                self.assertListEqual(
                    list(new),
                    list(old)
                )

        os.remove(new_path)

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

class TestXYZ(unittest.TestCase):
    def test_writefile(self):
        read_path = "test/static/test_peptide.xyz"
        path = "test/static/test_peptide.inp"
        new_path = "test/static/test_peptide_copy.inp"

        file = cctk.XYZFile.read_file(read_path)
        header = "! aug-cc-pVTZ aug-cc-pVTZ/C DLPNO-CCSD(T) TightSCF TightPNO MiniPrint\n%pal nproc 4 end\n%maxcore 4000\n%mdci\n    density none\nend"
        cctk.OrcaFile.write_molecule_to_file(new_path, file.molecule, header)

        with open(path) as old:
            with open(new_path) as new:
                self.assertListEqual(
                    list(new),
                    list(old)
                )

        os.remove(new_path)

        ensemble = cctk.ConformationalEnsemble()
        ensemble.add_molecule(file.molecule)

        orca_file = cctk.OrcaFile(ensemble=ensemble, header=header)
        orca_file.write_file(new_path)

        with open(path) as old:
            with open(new_path) as new:
                self.assertListEqual(
                    list(new),
                    list(old)
                )

        os.remove(new_path)

class TestPDB(unittest.TestCase):
    def test_write_traj(self):
        path = "test/static/gaussian_file.out"
        file = cctk.GaussianFile.read_file(path)
        mols = file.ensemble
        self.assertTrue(isinstance(mols, cctk.ConformationalEnsemble))

        old_path = "test/static/traj.pdb"
        new_path = "test/static/new_traj.pdb"

        cctk.PDBFile.write_ensemble_to_trajectory(new_path, mols)

        with open(old_path) as old:
            with open(new_path) as new:
                self.assertListEqual(
                    list(new),
                    list(old)
                )

        os.remove(new_path)

if __name__ == '__main__':
    unittest.main()
