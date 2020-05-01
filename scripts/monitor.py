import sys, re, glob, cctk
import numpy as np
import pandas as pd

from tabulate import tabulate
from asciichartpy import plot

#### This is a script to monitor the output of Gaussian files. 
#### In contrast to ``analyze.py``, this script analyzes only one file! 

#### Usage: ``python monitor.py path/to/output.out``

#### Corin Wagen and Eugene Kwan, 2019

filename = sys.argv[1]
print(f"\n\033[3mreading {filename}:\033[0m")

output_file = cctk.GaussianFile.read_file(filename, extended_opt_info=True)
if isinstance(output_file, list):
    output_file = output_file[-1]
print(f"{output_file.successful_terminations} successful terminations")
print(f"{output_file.num_imaginaries()} imaginary frequencies")
if output_file.num_imaginaries():
    freqs = [f"{f:.1f} cm-1" for f in output_file.imaginaries()]
    for f in freqs:
        print(f"\t{f}")

print("\n\033[3manalysis:\033[0m")
print(f"{len(output_file.ensemble)} iterations completed")
property_names = ["scf_iters", "energy", "quasiharmonic_gibbs_free_energy", "rms_force", "rms_displacement", "rms_gradient"]

properties = output_file.ensemble[:,property_names]
if len(output_file.ensemble) == 1:
    properties = [properties]
df = pd.DataFrame(properties, columns=property_names).fillna(0)
df["rel_energy"] = (df.energy - df.energy.min()) * 627.509469

print("\n\033[1mENERGY (kcal/mol):\033[0m")
print(plot(df["rel_energy"], {"height": 12, "format": " {:8.2f} "}))

print("\n\033[1mRMS GRADIENT:\033[0m")
print(plot(df["rms_gradient"], {"height": 12, "format": " {:8.6f} "}))

print("\n\033[1mRMS FORCE:\033[0m")
print(plot(df["rms_force"], {"height": 12, "format": " {:8.6f} "}))

print("\n\033[1mRMS DISPLACEMENT:\033[0m")
print(plot(df["rms_displacement"], {"height": 12, "format": " {:8.4f} "}))
