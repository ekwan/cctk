import unittest, sys, os, io, copy
import numpy as np
import cctk

class TestXYZ(unittest.TestCase):
    def test_readfile(self):
        path = "static/test_peptide.xyz"
        file = cctk.XYZFile.read_file(path)
        self.assertEqual(file.title, "peptide example")

        mol = file.molecule
        self.assertTrue(isinstance(mol, cctk.Molecule))
        self.assertEqual(mol.num_atoms(), 31)
        self.assertTrue(mol.check_for_conflicts())

    def test_writefile(self):
        path = "static/test_peptide.xyz"
        new_path = "static/test_peptide_copy.xyz"

        file = cctk.XYZFile.read_file(path)
        file.write_file(new_path)

        with open(path) as old:
            with open(new_path) as new:
                self.assertListEqual(
                    list(new),
                    list(old)
                )

        os.remove(new_path)

class TestGaussian(unittest.TestCase):
    def test_read_gjf_file(self):
        path = "static/gaussian_file.gjf"
        file = cctk.GaussianFile.read_file(path)
        self.assertEqual(file.header, "#p opt freq=noraman m062x/6-31g(d) scrf=(smd,solvent=diethylether)")
        self.assertEqual(file.title, "title")
        self.assertEqual(file.footer, None)

        mol = file.get_molecule()
        self.assertTrue(isinstance(mol, cctk.Molecule))
        self.assertEqual(mol.num_atoms(), 31)
        self.assertEqual(mol.charge, 0)
        self.assertEqual(mol.multiplicity, 1)

    def test_read_out_file(self):
        path = "static/gaussian_file.out"
        file = cctk.GaussianFile.read_file(path)
        self.assertEqual(file.header, "#p opt freq=noraman m062x/6-31g(d) scrf=(smd,solvent=diethylether)")
        self.assertEqual(file.title, "title")
        self.assertEqual(file.footer, None)
        self.assertTrue(isinstance(file.molecules, cctk.ConformationalEnsemble))

        mol = file.get_molecule()
        self.assertTrue(isinstance(mol, cctk.Molecule))
        self.assertEqual(mol.num_atoms(), 31)
        self.assertEqual(mol.charge, 0)
        self.assertEqual(mol.multiplicity, 1)

        old_path = "static/gaussian_file.gjf"
        new_path = "static/new_gjf.gjf"

        file.write_file(new_path, molecule=2)

        with open(old_path) as old:
            with open(new_path) as new:
                self.assertListEqual(
                    list(new),
                    list(old)
                )

        os.remove(new_path)

class TestMOL2(unittest.TestCase):
    def test_read(self):
        path = "static/dodecane.mol2"
        file = cctk.MOL2File.read_file(path)
        self.assertTrue(isinstance(file, cctk.MOL2File))

        ensemble = file.molecules
        self.assertTrue(isinstance(ensemble, cctk.ConformationalEnsemble))
        self.assertEqual(len(ensemble.molecules), 1)

        mol = ensemble.molecules[0]
        self.assertTrue(isinstance(mol, cctk.Molecule))
        self.assertEqual(len(mol.atomic_numbers), 38)
        self.assertEqual(len(mol.geometry), 38)
        self.assertEqual(mol.get_bond_order(1,2), 1)

    def test_bulk_read(self):
        path = "static/dodecane-csearch.mol2"
        file = cctk.MOL2File.read_file(path)
        self.assertTrue(isinstance(file, cctk.MOL2File))

        ensemble = file.molecules
        self.assertTrue(isinstance(ensemble, cctk.ConformationalEnsemble))
        self.assertEqual(len(ensemble.molecules), 597)
        for mol in ensemble.molecules:
            self.assertEqual(len(mol.atomic_numbers), 38)
            self.assertEqual(len(mol.geometry), 38)
            self.assertEqual(mol.get_bond_order(1,2), 1)

class TestMAE(unittest.TestCase):
    def test_read(self):
        path = "static/dodecane_csearch-out.mae"
        (file, pnames, pvals) = cctk.MAEFile.read_file(path)
        self.assertTrue(isinstance(file, cctk.MAEFile))
        self.assertEqual(len(pnames), 597)
        self.assertEqual(len(pvals), 597)

        ensemble = file.molecules
        self.assertTrue(isinstance(ensemble, cctk.ConformationalEnsemble))
        self.assertEqual(len(ensemble.molecules), 597)
        for mol in ensemble.molecules:
            self.assertEqual(len(mol.atomic_numbers), 38)
            self.assertEqual(len(mol.geometry), 38)
            self.assertEqual(mol.get_bond_order(1,2), 1)

class TestGroup(unittest.TestCase):
    def test_group_add(self):
        path = "static/acetaldehyde.out"
        old_path = "static/14-butanedione.gjf"
        new_path = "static/new_14-butanedione.gjf"

        file = cctk.GaussianFile.read_file(path)
        group = cctk.Group.new_from_molecule(attach_to=6, molecule=file.get_molecule())
        new_mol = cctk.Group.add_group_to_molecule(file.get_molecule(), group, 5)
        file.write_file("static/new_14-butanedione.gjf", molecule=new_mol)

        with open(old_path) as old:
            with open(new_path) as new:
                self.assertListEqual(
                    list(new),
                    list(old)
                )

        os.remove(new_path)

class TestXYZ(unittest.TestCase):
    def test_writefile(self):
        read_path = "static/test_peptide.xyz"
        path = "static/test_peptide.inp"
        new_path = "static/test_peptide_copy.inp"

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

        orca_file = cctk.OrcaFile(molecules=ensemble, header=header)
        orca_file.write_file(new_path)

        with open(path) as old:
            with open(new_path) as new:
                self.assertListEqual(
                    list(new),
                    list(old)
                )

        os.remove(new_path)

class TestEnsemble(unittest.TestCase):
    def generate_test_ensemble(self):
        path = "static/test_peptide.xyz"
        file = cctk.XYZFile.read_file(path)
        mol = file.molecule

        e1 = np.array([1, 0, 0])
        e2 = np.array([0, 1, 0])
        e3 = np.array([0, 0, 1])

        ensemble = cctk.ConformationalEnsemble()
        ensemble.add_molecule(mol)
        self.assertEqual(len(ensemble.molecules), 1)

        mol_rot = copy.deepcopy(ensemble.molecules[0]).rotate_molecule(e1, 90)
        ensemble.add_molecule(mol_rot)
        self.assertEqual(len(ensemble.molecules), 2)

        mol_trans = copy.deepcopy(ensemble.molecules[0]).translate_molecule(e2)
        ensemble.add_molecule(mol_trans)
        self.assertEqual(len(ensemble.molecules), 3)

        mol_trans_rot = copy.deepcopy(ensemble.molecules[1].translate_molecule(e2))
        ensemble.add_molecule(mol_trans_rot)
        self.assertEqual(len(ensemble.molecules), 4)

        mol_rot_trans = copy.deepcopy(ensemble.molecules[2].rotate_molecule(e1, 90))
        ensemble.add_molecule(mol_rot_trans)
        self.assertEqual(len(ensemble.molecules), 5)

        ensemble.add_molecule(copy.deepcopy(ensemble.molecules[4].rotate_molecule(e3, 20)))
        ensemble.add_molecule(copy.deepcopy(ensemble.molecules[5].rotate_molecule(e3, 20)))
        ensemble.add_molecule(copy.deepcopy(ensemble.molecules[6].rotate_molecule(e3, 20)))
        ensemble.add_molecule(copy.deepcopy(ensemble.molecules[7].rotate_molecule(e3, 20)))
        ensemble.add_molecule(copy.deepcopy(ensemble.molecules[4].rotate_molecule(e3, -20)))
        ensemble.add_molecule(copy.deepcopy(ensemble.molecules[5].rotate_molecule(e3, -20)))
        ensemble.add_molecule(copy.deepcopy(ensemble.molecules[6].rotate_molecule(e3, -20)))
        ensemble.add_molecule(copy.deepcopy(ensemble.molecules[7].rotate_molecule(e3, -20)))
        return ensemble

    def test_align(self):
        #### since all the molecules are identical, every way we do this should be totally fine
        ensemble = self.generate_test_ensemble()
        ensemble = ensemble.align()
        template = ensemble.molecules[0].geometry
        for molecule in ensemble.molecules:
            for i in range(1,len(template)+1):
                self.assertTrue(cctk.helper_functions.compute_distance_between(molecule.geometry[i],template[i]) < 0.0001)

        ensemble2 = self.generate_test_ensemble()
        ensemble2 = ensemble2.align(atoms="heavy")
        template = ensemble2.molecules[0].geometry
        for molecule in ensemble2.molecules:
            for i in range(1,len(template)+1):
                self.assertTrue(cctk.helper_functions.compute_distance_between(molecule.geometry[i],template[i]) < 0.0001)

        ensemble3 = self.generate_test_ensemble()
        ensemble3 = ensemble3.align(atoms=[13, 4, 27, 6, 9, 14])
        template = ensemble3.molecules[0].geometry
        for molecule in ensemble3.molecules:
            for i in range(1,len(template)+1):
                self.assertTrue(cctk.helper_functions.compute_distance_between(molecule.geometry[i],template[i]) < 0.0001)


if __name__ == '__main__':
    unittest.main()
