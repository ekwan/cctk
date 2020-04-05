import unittest, sys, os, io, copy
import numpy as np
import cctk

class TestGroup(unittest.TestCase):
    def test_group_add(self):
        path = "test/static/acetaldehyde.out"
        old_path = "test/static/14-butanedione.gjf"
        new_path = "test/static/new_14-butanedione.gjf"

        file = cctk.GaussianFile.read_file(path)
        group = cctk.Group.new_from_molecule(attach_to=6, molecule=file.get_molecule())
        new_mol = cctk.Group.add_group_to_molecule(file.get_molecule(), group, 5)
        file.write_file("test/static/new_14-butanedione.gjf", molecule=new_mol)

        with open(old_path) as old:
            with open(new_path) as new:
                self.assertListEqual(
                    list(new),
                    list(old)
                )

        os.remove(new_path)

if __name__ == '__main__':
    unittest.main()
