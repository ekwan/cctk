import unittest, sys, os, io, copy
import numpy as np
import cctk
import glob as glob

# tests ensemble.molecules indexing
#
# python -m unittest test.test_ensemble2.TestEnsemble2
class TestEnsemble(unittest.TestCase):
    def build_test_ensemble(self):
        path = "test/static/phenylpropane*.out"
        conformational_ensemble = cctk.ConformationalEnsemble()
        for filename in sorted(glob.glob(path)):
            gaussian_file = cctk.GaussianFile.read_file(filename)
            ensemble = gaussian_file.ensemble
            molecule = ensemble.molecules[-1]
            properties_dict = ensemble.get_properties_dict(molecule)
            conformational_ensemble.add_molecule(molecule,properties_dict)
        return conformational_ensemble

    def test_ensemble(self):
        conformational_ensemble = self.build_test_ensemble()

        m1 = conformational_ensemble.molecules[0]
        self.assertEqual(conformational_ensemble[m1,"filename"], 'test/static/phenylpropane_1.out')
        l1 = conformational_ensemble.molecules[0:2]
        self.assertEqual(len(l1), 2)
        m1 = l1[0]
        m2 = l1[1]
        self.assertEqual(conformational_ensemble[m1,"filename"], 'test/static/phenylpropane_1.out')
        self.assertEqual(conformational_ensemble[m2,"filename"], 'test/static/phenylpropane_2.out')
        l2 = conformational_ensemble.molecules[[0,2,3]]
        l3 = conformational_ensemble[l2,"filename"]
        self.assertListEqual(l3, ['test/static/phenylpropane_1.out', 'test/static/phenylpropane_3.out', 'test/static/phenylpropane_4.out'])
        l4 = conformational_ensemble.molecules[0:4:2]
        self.assertListEqual(conformational_ensemble[l4,"filename"], ['test/static/phenylpropane_1.out', 'test/static/phenylpropane_3.out'])
        m3 = conformational_ensemble.molecules[-1]
        self.assertEqual(conformational_ensemble[m3,"filename"], 'test/static/phenylpropane_6.out')
        with self.assertRaises(AssertionError):
            m4 = conformational_ensemble.molecules[-10]
        with self.assertRaises(AssertionError):
            m4 = conformational_ensemble.molecules[10]
        with self.assertRaises(AssertionError):
            m4 = conformational_ensemble.molecules[[1,7]]
        with self.assertRaises(ValueError):
            m4 = conformational_ensemble.molecules["abc"]

    def test_ensemble_indexing(self):
        path = "test/static/gaussian_file.out"
        file = cctk.GaussianFile.read_file(path)
        mols = file.ensemble
        self.assertTrue(isinstance(mols, cctk.ConformationalEnsemble))

        self.assertEqual(len(mols), 3)
        self.assertTrue(isinstance(mols[0], cctk.ConformationalEnsemble))

        self.assertListEqual(mols.get_property(None, "energy"), [-1159.56782625, -1159.56782622, -1159.56782622])
        self.assertListEqual(mols.get_property(None, "enthalpy"), [None, None, -1159.314817])

        self.assertListEqual(mols[:,"energy"], [-1159.56782625, -1159.56782622, -1159.56782622])
        self.assertListEqual(mols[:,"enthalpy"], [None, None, -1159.314817])

        self.assertEqual(mols[-1,"energy"], -1159.56782622)
        self.assertEqual(mols[-1,"enthalpy"], -1159.314817)
        self.assertEqual(mols[2, "energy"], -1159.56782622)
        self.assertEqual(mols[2, "enthalpy"], -1159.314817)

        mols[2, "potato"] = "russet"
        self.assertEqual(mols[2, "potato"], "russet")
        self.assertEqual(mols[-1, "potato"], "russet")
        mols[:, "oil_type"] = "grapeseed"
        self.assertEqual(mols[:, "oil_type"], ["grapeseed"] * 3) # nut allergies are no joke
        mols[1, ["colonel", "condiment"]] = "mustard"
        self.assertEqual(mols[1, "condiment"], "mustard")
        self.assertEqual(mols[1, "colonel"], "mustard") # cf. Clue (1985)

        self.assertListEqual(list(mols.keys()), list(mols._items.keys()))

        mols[1,"energy"] = 300
        self.assertListEqual(mols[0:10,"energy"], [-1159.56782625, 300, -1159.56782622])
        mols[1,"energy"] = 200
        self.assertListEqual(mols[0:10,"energy"], [-1159.56782625, 200, -1159.56782622])
        mols[[1,2],"energy"] = [200, 201]
        self.assertListEqual(mols[0:10,"energy"], [-1159.56782625, 200, 201])
        mols[1:3,"energy"] = [203, 204]
        self.assertListEqual(mols[0:10,"energy"], [-1159.56782625, 203, 204])
        mols[[1,2],"energy"] = [100, 101]
        self.assertListEqual(mols[0:10,"energy"], [-1159.56782625, 100, 101])

        new_ensemble = mols[[1,2]]
        self.assertEqual(len(new_ensemble), 2)

        index_ndarray = np.array([1,2])
        new_ensemble = mols[index_ndarray]
        self.assertEqual(len(new_ensemble), 2)

        for (m, p) in mols:
            self.assertTrue(isinstance(m, cctk.Molecule))
            self.assertTrue(isinstance(p, dict))

    def test_ensemble_indexing2(self):
        path = "test/static/gaussian_file.out"
        file = cctk.GaussianFile.read_file(path)
        mols = file.ensemble
        self.assertTrue(isinstance(mols, cctk.ConformationalEnsemble))

        ensemble = cctk.ConformationalEnsemble()
        for i, molecule in enumerate(mols.molecules):
            ensemble.add_molecule(molecule)
            ensemble[molecule,"test_property"]=i

        self.assertTrue(len(ensemble)==3)
        self.assertListEqual(list(ensemble[:,"test_property"]), [0, 1, 2])

        ensemble = cctk.Ensemble()
        for i,molecule in enumerate(mols.molecules):
            ensemble.add_molecule(molecule)
            ensemble[molecule,"test_property"]=i

        self.assertTrue(len(ensemble)==3)
        self.assertListEqual(list(ensemble[:,"test_property"]), [0, 1, 2])

    def test_ensemble_properties(self):
        filename = "test/static/gaussian_file.out"
        gaussian_file = cctk.GaussianFile.read_file(filename)
        ensemble = gaussian_file.ensemble
        self.assertTrue(isinstance(ensemble, cctk.ConformationalEnsemble))

        self.assertEqual(len(ensemble), 3)
        self.assertTrue(isinstance(ensemble[0], cctk.ConformationalEnsemble))

        self.assertEqual(ensemble[0,"energy"], -1159.56782625)
        self.assertListEqual(ensemble[:,"energy"], [-1159.56782625, -1159.56782622, -1159.56782622])
        self.assertListEqual(ensemble[:,"enthalpy"], [None, None, -1159.314817])

    def test_sort(self):
        conformational_ensemble = self.build_test_ensemble()

        original_order = conformational_ensemble[:,"energy"]
        self.assertListEqual(original_order,[0.0140132996483, 0.0163679933924, 0.0213666533731, 0.0180903133947, 0.0547890926923, 0.0182782865186])
        sorted_ensemble = conformational_ensemble.sort_by("energy", ascending=False)
        sorted_order = sorted_ensemble[:,"energy"]
        self.assertListEqual(sorted_order,[0.0547890926923, 0.0213666533731, 0.0182782865186, 0.0180903133947, 0.0163679933924, 0.0140132996483])
        sorted_ensemble[2,"energy"]=None
        with self.assertRaises(ValueError):
            sorted_ensemble = conformational_ensemble.sort_by("energy", ascending=False)
        sorted_ensemble[2,"energy"]=0.0182782865186
        lowest_energy_molecules = conformational_ensemble.lowest_molecules("energy",2)
        self.assertEqual(len(lowest_energy_molecules), 2)
        energy0 = sorted_ensemble.get_property(lowest_energy_molecules[0], "energy")
        self.assertEqual(energy0, 0.0140132996483)
        energy1 = sorted_ensemble.get_property(lowest_energy_molecules[1], "energy")
        self.assertEqual(energy1, 0.0163679933924)
        lowest_molecule = conformational_ensemble.lowest_molecules("energy",1)
        self.assertTrue(isinstance(lowest_molecule, cctk.Molecule))
        energy0 = sorted_ensemble.get_property(lowest_molecule, "energy")
        self.assertEqual(energy0, 0.0140132996483)

    def test_boltzmann_weighting(self):
        conformational_ensemble = self.build_test_ensemble()

        values, weights = conformational_ensemble.boltzmann_average("energy", energies=[1.36,0,1000,1000,1000,1000], energy_unit="kcal_mol", return_weights=True)
        self.assertTrue((weights[0]/weights[1] - 0.1 < 0.01))
        self.assertTrue(values - 0.016152 < 0.0001)

        ce2 = cctk.ConformationalEnsemble()
        for filename in glob.glob("test/static/pentane*.out"):
            gaussian_file = cctk.GaussianFile.read_file(filename)
            ensemble = gaussian_file.ensemble
            molecule = ensemble.molecules[-1]
            properties_dict = ensemble.get_properties_dict(molecule)
            ce2.add_molecule(molecule,properties_dict)
        enthalpy = ce2.boltzmann_average("enthalpy")
        self.assertTrue(enthalpy - .10722 < 0.0001)

if __name__ == '__main__':
    unittest.main()
