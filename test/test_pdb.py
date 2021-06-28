import unittest, sys, os
import cctk

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
