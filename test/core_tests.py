import unittest, sys, os, io
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

class TestMolecule(unittest.TestCase):
    def load_molecule(self, path="static/test_peptide.xyz"):
        return cctk.XYZFile.read_file(path).molecule

    def test_distance(self):
        mol = self.load_molecule()
        self.assertTrue(isinstance(mol, cctk.Molecule))

        self.assertEqual(int(round(mol.get_distance(1,2)*10)), 10)
        self.assertEqual(int(round(mol.get_distance(1,3)*10)), 14)
        self.assertEqual(int(round(mol.get_distance(1,9)*10)), 38)

        mol.set_distance(1, 2, 2.00)

        self.assertEqual(int(round(mol.get_distance(1,2)*10)), 20)
        self.assertEqual(int(round(mol.get_distance(1,3)*10)), 14)
        self.assertEqual(int(round(mol.get_distance(1,9)*10)), 38)

        self.assertTrue(mol.check_for_conflicts())
        mol.set_distance(1, 2, 0.01)
        self.assertRaises(ValueError, mol.check_for_conflicts)

    def test_angle(self):
        mol = self.load_molecule()
        self.assertTrue(isinstance(mol, cctk.Molecule))

        self.assertEqual(int(round(mol.get_angle(1,3,5))), 111)
        self.assertEqual(int(round(mol.get_angle(3,5,7))), 110)
        self.assertEqual(int(round(mol.get_angle(5,7,9))), 64)

        mol.set_angle(1, 3, 5, 120)

        self.assertEqual(int(round(mol.get_angle(1,3,5))), 120)

    def test_dihedral(self):
        mol = self.load_molecule()
        self.assertTrue(isinstance(mol, cctk.Molecule))

        self.assertEqual(int(round(mol.get_dihedral(1,3,5,7))), 60)
        self.assertEqual(int(round(mol.get_dihedral(16,14,17,18))), 11)
        self.assertEqual(int(round(mol.get_dihedral(31,28,1,2))), 148)

        mol.set_dihedral(1, 3, 5, 7, 120)

        self.assertEqual(int(round(mol.get_dihedral(1,3,5,7))), 120)

    def test_fragment(self):
        mol = self.load_molecule()
        mol.assign_connectivity()
        (frag1, frag2) = mol._get_bond_fragments(3, 5)
        self.assertEqual(len(frag1), 27)
        self.assertEqual(len(frag2), 4)

        self.assertEqual(len(mol._get_fragment_containing(5)), 31)
        mol.remove_bond(3,5)
        self.assertEqual(len(mol._get_fragment_containing(5)), 4)
        self.assertFalse(mol.are_connected(3,5))
        mol.add_bond(3,5)
        self.assertEqual(len(mol._get_fragment_containing(5)), 31)
        self.assertTrue(mol.are_connected(3,5))

    def test_add_atoms(self):
        mol = cctk.Molecule(np.array([2], dtype=np.int8), [[0, 0, 0]])
        self.assertEqual(mol.num_atoms(), 1)

        mol.add_atom("He", [1, 0, 0])
        self.assertListEqual(list(mol.atomic_numbers), [2, 2])
        self.assertEqual(mol.num_atoms(), 2)

        mol.add_atom("Ar", [3, 0, 0])
        self.assertEqual(mol.num_atoms(), 3)

        mol.add_atom_at_centroid("He", [2, 3])
        self.assertEqual(mol.num_atoms(), 4)
        self.assertListEqual(list(mol.get_vector(4)), [2, 0, 0])

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

if __name__ == '__main__':
    unittest.main()
