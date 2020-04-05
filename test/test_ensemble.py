import unittest, sys, os, io, copy
import numpy as np
import cctk

class TestEnsemble(unittest.TestCase):
    def generate_test_ensemble(self):
        path = "test/static/test_peptide.xyz"
        file = cctk.XYZFile.read_file(path)
        mol = file.molecule

        e1 = np.array([1, 0, 0])
        e2 = np.array([0, 1, 0])
        e3 = np.array([0, 0, 1])

        ensemble = cctk.ConformationalEnsemble()
        ensemble.add_molecule(mol)
        self.assertEqual(len(ensemble), 1)

        mol_rot = copy.deepcopy(ensemble[0]).rotate_molecule(e1, 90)
        ensemble.add_molecule(mol_rot)
        self.assertEqual(len(ensemble), 2)

        mol_trans = copy.deepcopy(ensemble[0]).translate_molecule(e2)
        ensemble.add_molecule(mol_trans)
        self.assertEqual(len(ensemble), 3)

        mol_trans_rot = copy.deepcopy(ensemble[1].translate_molecule(e2))
        ensemble.add_molecule(mol_trans_rot)
        self.assertEqual(len(ensemble), 4)

        mol_rot_trans = copy.deepcopy(ensemble[2].rotate_molecule(e1, 90))
        ensemble.add_molecule(mol_rot_trans)
        self.assertEqual(len(ensemble), 5)

        ensemble.add_molecule(copy.deepcopy(ensemble[4].rotate_molecule(e3, 20)))
        ensemble.add_molecule(copy.deepcopy(ensemble[5].rotate_molecule(e3, 20)))
        ensemble.add_molecule(copy.deepcopy(ensemble[6].rotate_molecule(e3, 20)))
        ensemble.add_molecule(copy.deepcopy(ensemble[7].rotate_molecule(e3, 20)))
        ensemble.add_molecule(copy.deepcopy(ensemble[4].rotate_molecule(e3, -20)))
        ensemble.add_molecule(copy.deepcopy(ensemble[5].rotate_molecule(e3, -20)))
        ensemble.add_molecule(copy.deepcopy(ensemble[6].rotate_molecule(e3, -20)))
        ensemble.add_molecule(copy.deepcopy(ensemble[7].rotate_molecule(e3, -20)))
        return ensemble

    def test_align(self):
        #### since all the molecules are identical, every way we do this should be totally fine
        ensemble = self.generate_test_ensemble()
        ensemble = ensemble.align()
        template = ensemble[0].geometry
        for molecule in ensemble.molecules():
            for i in range(1,len(template)+1):
                self.assertTrue(cctk.helper_functions.compute_distance_between(molecule.geometry[i],template[i]) < 0.0001)

        ensemble2 = self.generate_test_ensemble()
        ensemble2 = ensemble2.align(comparison_atoms="heavy")
        template = ensemble2[0].geometry
        for molecule in ensemble2.molecules():
            for i in range(1,len(template)+1):
                self.assertTrue(cctk.helper_functions.compute_distance_between(molecule.geometry[i],template[i]) < 0.0001)

        ensemble3 = self.generate_test_ensemble()
        ensemble3 = ensemble3.align(comparison_atoms=[13, 4, 27, 6, 9, 14])
        template = ensemble3[0].geometry
        for molecule in ensemble3.molecules():
            for i in range(1,len(template)+1):
                self.assertTrue(cctk.helper_functions.compute_distance_between(molecule.geometry[i],template[i]) < 0.0001)

if __name__ == '__main__':
    unittest.main()
