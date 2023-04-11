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
        self.assertTrue(isinstance(file.get_molecule(), cctk.Molecule))

        header = "! aug-cc-pVTZ aug-cc-pVTZ/C DLPNO-CCSD(T) TightSCF TightPNO"
        variables = {"maxcore": 4000}
        blocks = {"pal": ["nproc 4"], "mdci": ["density none"]}

        cctk.OrcaFile.write_molecule_to_file(new_path, file.get_molecule(), header, variables, blocks)

        with open(path) as old:
            with open(new_path) as new:
                self.assertListEqual(
                    list(new),
                    list(old)
                )

        os.remove(new_path)

        ensemble = cctk.ConformationalEnsemble()
        ensemble.add_molecule(file.get_molecule())

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
        path = "test/static/H2O_dlpno_ccsdt.out"
        file = cctk.OrcaFile.read_file(path)
        self.assertEqual(file.successful_terminations, 1)
        self.assertEqual(file.elapsed_time, 18.471)
        self.assertEqual(file.header, "! cc-pVTZ cc-pVTZ/C DLPNO-CCSD(T) TightSCF TightPNO")
        self.assertEqual(file.variables["maxcore"], "1000")
        self.assertListEqual(file.blocks["mdci"], ["density none"])

        mol = file.get_molecule()
        self.assertTrue(isinstance(mol, cctk.Molecule))
        self.assertEqual(mol.num_atoms(), 3)
        self.assertEqual(file.ensemble[mol,"energy"], -76.330947653965)

        path = "test/static/AcOH_orca.out"
        file = cctk.OrcaFile.read_file(path)

        mol = file.get_molecule()
        self.assertTrue(isinstance(mol, cctk.Molecule))
        self.assertEqual(mol.num_atoms(), 8)
        self.assertEqual(file.ensemble[mol,"energy"], -229.12132242363)
        self.assertEqual(file.ensemble[mol, "dipole_moment"], 1.76241)

        self.assertEqual(file.ensemble[mol, "mulliken_charges"][1], 0.329311)
        self.assertEqual(file.ensemble[mol, "lowdin_charges"][1], -0.539274)

        self.assertEqual(file.ensemble[mol, "temperature"], 298.15)
        self.assertEqual(file.ensemble[mol, "enthalpy"], -229.05330337)
        self.assertEqual(file.ensemble[mol, "gibbs_free_energy"], -229.08534132)
        self.assertListEqual(list(file.ensemble[mol, 'frequencies'][:3]), [129.95, 432.62, 559.79])

    def test_nmr(self):
        path = "test/static/ibuprofen_nmr_orca.out"
        file = cctk.OrcaFile.read_file(path)

        molecule = file.get_molecule()
        properties_dict = file.ensemble.get_properties_dict(molecule)
        energy = properties_dict["energy"]
        self.assertTrue(abs(energy + 656.306067336866) < 1e8)
        self.assertTrue(abs(file.ensemble[-1, "energy"] + 656.306067336866) < 1e8)
        shieldings = properties_dict["isotropic_shielding"]
        self.assertListEqual(list(shieldings[:6]), [55.307, 68.003, 63.738, 51.446, 65.325])

    def test_freq(self):
        # this is a transition state search
        path = "test/static/orca_OptTs.out"
        file = cctk.OrcaFile.read_file(path)
        freqs = file.ensemble[-1, 'frequencies']
        self.assertListEqual(list(freqs[:3]), [-2742.92, -25.93, 386.57])

        # this is a transition state search
        # that recalculates the Hessian every 2 steps
        path = "test/static/orca_OptTs_RecalcHess.out"
        file = cctk.OrcaFile.read_file(path)
        freqs = file.ensemble[-1, 'frequencies']
        self.assertListEqual(list(freqs[:3]), [-2742.91, -25.92, 386.56])

    def test_compound_job(self):
        # this is a compound job with opt/freq then sp.
        path = "test/static/orca_gemfi_alfa_minima_1.out"
        files = cctk.OrcaFile.read_file(path)

        # test the properties of the opt/freq job
        file = files[0]
        mol = file.get_molecule()
        self.assertEqual(file.ensemble[mol,"energy"], -1493.894058726923)
        self.assertEqual(file.ensemble[mol, "dipole_moment"], 28.00267)
        self.assertEqual(file.ensemble[mol, "mulliken_charges"][1], -0.840055)
        self.assertEqual(file.ensemble[mol, "lowdin_charges"][1], -0.024951)
        self.assertEqual(file.ensemble[mol, "temperature"], 298.15)
        self.assertEqual(file.ensemble[mol, "enthalpy"], -1493.38539502)
        self.assertListEqual(list(file.ensemble[mol, 'frequencies'][:3]), [-55.15, 16.13, 20.71])

        # test the properties of the subsequent single point calculation
        file = files[1]
        mol = file.get_molecule()
        self.assertEqual(file.ensemble[mol,"energy"], -1494.296011093052)
        self.assertEqual(file.ensemble[mol, "dipole_moment"], 28.10263)
        self.assertEqual(file.ensemble[mol, "mulliken_charges"][1], -0.630256)
        self.assertEqual(file.ensemble[mol, "lowdin_charges"][1], -0.044225)

        # this is a compound job with ScanTs/freq followed by sp.
        path = "test/static/orca_gemfi_beta_dTS_1.out"
        files = cctk.OrcaFile.read_file(path)

        # test the properties of the ScanTs/freq job
        file = files[0]
        mol = file.get_molecule()
        self.assertEqual(file.ensemble[mol,"energy"], -1493.878871250123)
        # self.assertEqual(file.ensemble[mol, "dipole_moment"], 36.87591) # if we update to support parsing scan jobs
        self.assertEqual(file.ensemble[mol, "mulliken_charges"][1], 0.002160)
        self.assertEqual(file.ensemble[mol, "enthalpy"], -1493.37176571)
        self.assertListEqual(list(file.ensemble[mol, 'frequencies'][:3]), [-213.88, 13.56, 19.58])

        file = files[1]
        mol = file.get_molecule()
        self.assertEqual(file.ensemble[mol,"energy"], -1494.275051149283)
        self.assertEqual(file.ensemble[mol, "dipole_moment"], 37.12405)
        self.assertEqual(file.ensemble[mol, "mulliken_charges"][1], -0.211305)
        self.assertEqual(file.ensemble[mol, "lowdin_charges"][1], -0.046969)








        








