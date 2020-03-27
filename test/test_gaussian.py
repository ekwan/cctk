import unittest, sys, os, io, copy
import numpy as np
import cctk

class TestGaussian(unittest.TestCase):
    def test_read_gjf_file(self):
        path = "test/static/gaussian_file.gjf"
        file = cctk.GaussianFile.read_file(path)
        self.assertEqual(file.route_card, "#p opt freq=noraman m062x/6-31g(d) scrf=(smd,solvent=diethylether)")
        self.assertListEqual(file.job_types, [cctk.JobType.OPT, cctk.JobType.FREQ])
        self.assertDictEqual(file.link0, {"mem": "1GB", "chk": "test.chk"})
        self.assertEqual(file.title, "title")
        self.assertEqual(file.footer, None)

        mol = file.get_molecule()
        self.assertTrue(isinstance(mol, cctk.Molecule))
        self.assertEqual(mol.num_atoms(), 31)
        self.assertEqual(mol.charge, 0)
        self.assertEqual(mol.multiplicity, 1)

    def test_read_out_file(self):
        path = "test/static/gaussian_file.out"
        file = cctk.GaussianFile.read_file(path)
        self.assertEqual(file.route_card, "#p opt freq=noraman m062x/6-31g(d) scrf=(smd,solvent=diethylether)")
        self.assertDictEqual(file.link0, {"mem": "32GB",  "nprocshared": "16"})
        self.assertListEqual(file.job_types, [cctk.JobType.OPT, cctk.JobType.FREQ])
        self.assertEqual(file.title, "title")
        self.assertEqual(file.footer, None)
        self.assertTrue(isinstance(file.molecules, cctk.ConformationalEnsemble))

        mol = file.get_molecule()
        self.assertTrue(isinstance(mol, cctk.Molecule))
        self.assertEqual(mol.num_atoms(), 31)
        self.assertEqual(mol.charge, 0)
        self.assertEqual(mol.multiplicity, 1)

        old_path = "test/static/gaussian_file.gjf"
        new_path = "test/static/new_gjf.gjf"

        file.write_file(new_path, molecule=1, link0={"mem": "1GB", "chk": "test.chk"})

        with open(old_path) as old:
            with open(new_path) as new:
                self.assertListEqual(
                    list(new),
                    list(old)
                )

        os.remove(new_path)

    def test_link1_out_file(self):
        path = "test/static/ethane.out"
        f, lines = cctk.GaussianFile.read_file(path, return_lines=True)

        self.assertEqual(len(lines), 3)
        self.assertEqual(len(f), 3)
        self.assertTrue(all(isinstance(file, cctk.GaussianFile) for file in f))

        self.assertListEqual(f[0].job_types, [cctk.JobType.OPT, cctk.JobType.FREQ])
        self.assertListEqual(f[1].job_types, [cctk.JobType.NMR])
        self.assertListEqual(f[2].job_types, [cctk.JobType.NMR])

if __name__ == '__main__':
    unittest.main()
