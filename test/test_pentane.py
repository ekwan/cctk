import unittest
import numpy as np
import cctk
import glob as glob
import pandas as pd
from pandas import DataFrame

class TestPentane(unittest.TestCase):
    def test_penane(self):
        path = "test/static/pentane_conformation*.out"
        conformational_ensemble = cctk.ConformationalEnsemble()
        filenames = sorted(glob.glob(path))
        for filename in filenames:
            gaussian_file = cctk.GaussianFile.read_file(filename)
            ensemble = gaussian_file.ensemble
            molecule = ensemble.molecules[-1]
            properties_dict = ensemble.get_properties_dict(molecule)
            conformational_ensemble.add_molecule(molecule,properties_dict)
        property_names = ["filename", "energy"]
        conformational_energies = conformational_ensemble[:,property_names]
        df = DataFrame(conformational_energies, columns=property_names)
        df["rel_energy"] = (df.energy - df.energy.min()) * 627.509469
        #print(df)

if __name__ == '__main__':
    unittest.main()
