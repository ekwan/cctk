import sys, re, glob, cctk
import numpy as np
import pandas as pd

#### This is a script to monitor the output of Gaussian files. 
#### In contrast to ``monitor.py``, this script analyzes many files! 
#### If the file has not successfully achieved SCF convergence at least once, that file will not display any information. 

#### Usage: ``python analyze.py "path/to/output/*.out"``
#### NOTE: It's crucial to wrap the wildcard-containing path in quotes!

#### NOTE: This file will reject any file that contains the string "slurm."

#### Corin Wagen and Eugene Kwan, 2019

filenames = sys.argv[1]
results = cctk.Ensemble()

for filename in sorted(glob.glob(filenames, recursive=True)):
    print(f"reading file {filename}")
    if re.search("slurm", filename):
        continue

    output_file = cctk.GaussianFile.read_file(filename)
    if isinstance(output_file, list):
        output_file = output_file[-1]

    results.add_molecule(*list(output_file.ensemble.items())[-1])
    results[output_file.get_molecule(), "iterations"] = len(output_file.ensemble)
    results[output_file.get_molecule(), "success"] = output_file.success
    results[output_file.get_molecule(), "num_imag"] = output_file.num_imaginaries()

if len(results) == 0:
    print("no jobs to analyze!")
    exit()

property_names = ["filename", "iterations", "energy", "enthalpy", "gibbs_free_energy", "rms_force", "rms_displacement", "success", "num_imag"]
values = results[:, property_names]
if not isinstance(values[0], list):
    values = [values]

df = pd.DataFrame(values, columns=property_names)
df["rel_energy"] = (df.energy - df.energy.min()) * 627.509469
print(df)

