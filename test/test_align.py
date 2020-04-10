import unittest, sys, os, io, copy
import numpy as np
import cctk
import glob as glob

# python -m unittest test.test_align.TestAlign
class TestAlign(unittest.TestCase):
    def test_RMSD(self):
        path = "test/static/gaussian_file.out"
        gaussian_file = cctk.GaussianFile.read_file(path)
        conformational_ensemble = gaussian_file.molecules
        m1 = conformational_ensemble[0]
        m2 = conformational_ensemble[-1]
        RMSD = cctk.helper_functions.compute_RMSD(m1,m2)
        delta = abs(0.0006419131435567976 - RMSD)
        self.assertLess(delta, 0.0001)

    def test_align(self):
        path = "test/static/phenylpropane*.out"
        conformational_ensemble = cctk.ConformationalEnsemble()
        for filename in sorted(glob.glob(path)):
            gaussian_file = cctk.GaussianFile.read_file(filename)
            e = gaussian_file.molecules
            m = e[-1]
            p = e[m]
            conformational_ensemble.add_molecule(m,p)

        comparison_atoms = [1,2,3,4,5,6]
        # added np.int64 here to check that the to_geometry parameter will take any int
        aligned_ensemble, before_RMSD, after_RMSD = conformational_ensemble.align(to_geometry=np.int64(0), comparison_atoms=comparison_atoms, compute_RMSD=True)
        for before,after in zip(before_RMSD, after_RMSD):
            self.assertLess(after,0.0001)
        cctk.GaussianFile.write_ensemble_to_file("test/static/phenylpropane_aligned.gjf", aligned_ensemble, "#p")
        #for molecule,properties in aligned_ensemble.items():
        #    for k,v in properties.items():
        #        print(k,v)
        #    print
        #ensemble2 = aligned_ensemble.eliminate_redundant()

        #for i,molecule in enumerate(aligned_ensemble):
        #    filename = f"test/static/phenylpropane_{i+1}.gjf"
        #    route_card = "#"
        #    cctk.GaussianFile.write_molecule_to_file(filename, molecule, route_card)
        #    rmsd = cctk.helper_functions.compute_RMSD(m0,molecule,comparison_atoms)
        #    #print("after:",rmsd)

if __name__ == '__main__':
    unittest.main()
