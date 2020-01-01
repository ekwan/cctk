import unittest, sys, os, io
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

        mol = file.get_molecule()
        self.assertTrue(isinstance(mol, cctk.Molecule))
        self.assertEqual(mol.num_atoms(), 31)
        self.assertEqual(mol.charge, 0)
        self.assertEqual(mol.multiplicity, 1)

    def test_read_out_file(self):
        path = "static/gaussian_file.out"
        file = cctk.GaussianFile.read_file(path)

        mol = file.get_molecule()

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

if __name__ == '__main__':
    unittest.main()
