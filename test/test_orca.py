import unittest, sys, os, io, copy
import numpy as np
import cctk

if __name__ == '__main__':
    unittest.main()

class TestOrca(unittest.TestCase):
    def test_write(self):
        read_path = "test/static/test_peptide.xyz"
        path = "test/static/test_peptide.inp"
        new_path = "test/static/test_peptide_copy.inp"

        file = cctk.XYZFile.read_file(read_path)
        self.assertTrue(isinstance(file.molecule, cctk.Molecule))

        header = "! aug-cc-pVTZ aug-cc-pVTZ/C DLPNO-CCSD(T) TightSCF TightPNO MiniPrint"
        variables = {"maxcore": 4000}
        blocks = {"pal": ["nproc 4"], "mdci": ["density none"]}

        cctk.OrcaFile.write_molecule_to_file(new_path, file.molecule, header, variables, blocks)

        with open(path) as old:
            with open(new_path) as new:
                self.assertListEqual(
                    list(new),
                    list(old)
                )

        os.remove(new_path)

        ensemble = cctk.ConformationalEnsemble()
        ensemble.add_molecule(file.molecule)

        orca_file = cctk.OrcaFile(job_types=[cctk.OrcaJobType.SP], ensemble=ensemble, header=header, blocks=blocks, variables=variables)
        orca_file.write_file(new_path)

        with open(path) as old:
            with open(new_path) as new:
                self.assertListEqual(
                    list(new),
                    list(old)
                )

        os.remove(new_path)

    def test_read(self):
        path = "test/static/MsOH_ccsdt.out"
        file = cctk.OrcaFile.read_file(path)
        self.assertEqual(file.successful_terminations, 1)
        self.assertEqual(file.elapsed_time, 8575)
        self.assertEqual(file.header, "! aug-cc-pVQZ aug-cc-pVQZ/C DLPNO-CCSD(T) TightSCF TightPNO MiniPrint")
        self.assertEqual(file.variables["maxcore"], "50000")
        self.assertListEqual(file.blocks["mdci"], ["density none"])

        mol = file.get_molecule()
        self.assertTrue(isinstance(mol, cctk.Molecule))
        self.assertEqual(mol.num_atoms(), 9)
        self.assertEqual(file.ensemble[mol,"energy"], -663.663569902734)

        path = "test/static/AcOH_orca.out"
        file = cctk.OrcaFile.read_file(path)

        mol = file.get_molecule()
        self.assertTrue(isinstance(mol, cctk.Molecule))
        self.assertEqual(mol.num_atoms(), 8)
        self.assertEqual(file.ensemble[mol,"energy"], -229.12132242363)
        self.assertEqual(file.ensemble[mol, "dipole_moment"], 1.76241)

        self.assertEqual(file.ensemble[mol, "mulliken_charges"][1], 0.333851)
        self.assertEqual(file.ensemble[mol, "lowdin_charges"][1], -0.515118)

        self.assertEqual(file.ensemble[mol, "temperature"], 298.15)
        self.assertEqual(file.ensemble[mol, "enthalpy"], -229.05330337)
        self.assertEqual(file.ensemble[mol, "gibbs_free_energy"], -229.08534132)

