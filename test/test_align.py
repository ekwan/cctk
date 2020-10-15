import unittest, sys, os, io, copy
import numpy as np
import cctk
import glob as glob

# python -m unittest test.test_align.TestAlign
class TestAlign(unittest.TestCase):
    def test_RMSD(self):
        path = "test/static/gaussian_file.out"
        gaussian_file = cctk.GaussianFile.read_file(path)
        ensemble = gaussian_file.ensemble
        m1 = ensemble.molecules[0]
        m2 = ensemble.molecules[-1]
        RMSD = cctk.helper_functions.compute_RMSD(m1.geometry,m2.geometry)
        delta = abs(0.0006419131435567976 - RMSD)
        self.assertLess(delta, 0.0001)

    def test_align(self):
        path = "test/static/phenylpropane*.out"
        conformational_ensemble = cctk.ConformationalEnsemble()
        for filename in sorted(glob.glob(path)):
            gaussian_file = cctk.GaussianFile.read_file(filename)
            ensemble = gaussian_file.ensemble
            molecule = ensemble.molecules[-1]
            properties_dict = ensemble.get_properties_dict(molecule)
            conformational_ensemble.add_molecule(molecule,properties_dict)

        comparison_atoms = [1,2,3,4,5,6]
        # added np.int64 here to check that the to_geometry parameter will take any int
        aligned_ensemble, before_RMSD, after_RMSD = conformational_ensemble.align(to_geometry=np.int64(0), comparison_atoms=comparison_atoms, compute_RMSD=True)
        for before,after in zip(before_RMSD, after_RMSD):
            self.assertLess(after,0.0001)
        cctk.GaussianFile.write_ensemble_to_file("test/static/phenylpropane_aligned.gjf", aligned_ensemble, "#p")

        ensemble2, rmsds = aligned_ensemble.eliminate_redundant(RMSD_cutoff=0.5, comparison_atoms="heavy", return_RMSD=True)
        self.assertEqual(len(ensemble2), 3)

        cctk.GaussianFile.write_ensemble_to_file("test/static/phenylpropane_aligned2.gjf", ensemble2, "#p")
        ensemble3 = aligned_ensemble.eliminate_redundant(RMSD_cutoff=0.5, comparison_atoms=comparison_atoms)
        self.assertEqual(len(ensemble3), 1)

        cctk.GaussianFile.write_ensemble_to_file("test/static/phenylpropane_aligned3.gjf", ensemble3, "#p")
        cctk.MOL2File.write_ensemble_to_file("test/static/phenylpropane_aligned.mol2", aligned_ensemble)


if __name__ == '__main__':
    unittest.main()
