import unittest, sys, os, io, copy
import numpy as np
import cctk

class TestXYZ(unittest.TestCase):
    def test_readfile(self):
        path = "test/static/test_peptide.xyz"
        file = cctk.XYZFile.read_file(path)
        self.assertEqual(file.titles[0], "peptide example")

        mol = file.get_molecule()
        self.assertTrue(isinstance(mol, cctk.Molecule))
        self.assertEqual(mol.num_atoms(), 31)
        self.assertTrue(mol.check_for_conflicts())

        with self.assertWarns(DeprecationWarning):
            self.assertTrue(isinstance(file.molecule, cctk.Molecule))

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

    def test_traj(self):
        path = "test/static/methane_traj.xyz"
        file = cctk.XYZFile.read_trajectory(path)

        self.assertEqual(len(file.ensemble), 251)
        self.assertEqual(file.get_molecule().num_atoms(), 5)

    def test_ense(self):
        path = "test/static/methane_traj.xyz"
        file = cctk.XYZFile.read_ensemble(path)
        self.assertEqual(len(file.ensemble), 251)

        new_path = "test/static/methane_traj_new.xyz"
        cctk.XYZFile.write_ensemble_to_file(new_path, file.ensemble, title="sample title")
        os.remove(new_path)

if __name__ == '__main__':
    unittest.main()
