import unittest, sys, os, io, copy, math
import numpy as np
import networkx as nx
import cctk

class TestSolventExtraction(unittest.TestCase):
    def test_basic(self):
        mol = cctk.XYZFile.read_file("test/static/acetone_water.xyz").get_molecule()
        mol.assign_connectivity()
        self.assertFalse(mol.num_atoms() == 40)
        new_mol = mol.limit_solvent_shell(num_solvents=10)
        self.assertTrue(new_mol.num_atoms() == 40)

        path = "test/static/acetone_10waters.gjf"
        new_path = "test/static/new_acetone_10waters.gjf"

        with self.assertWarns(Warning):
            cctk.GaussianFile.write_molecule_to_file(
                new_path,
                new_mol,
                link0={"nprocshared": 4, "mem": "3GB"},
                route_card="#t b3lyp empiricaldispersion=gd3bj gen NMR pop=none int=finegrid nosymm",
                footer="@/n/jacobsen_lab/ekwan/solvent/basis/pcSseg-1.bas",
            )

        with open(path) as old:
            with open(new_path) as new:
                self.assertListEqual(
                    list(new),
                    list(old)
                )

        os.remove(new_path)

if __name__ == '__main__':
    unittest.main()
