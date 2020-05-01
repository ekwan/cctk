import sys, re, glob, cctk, argparse
import numpy as np
import pandas as pd

from tabulate import tabulate
from tqdm import tqdm

#### This is a script to monitor the output of Orca files. 
#### In contrast to ``monitor_orca.py``, this script analyzes many files! 
#### If the file has not successfully achieved SCF convergence at least once, that file will not display any information. 

#### Usage: ``python analyze_orca.py "path/to/output/*.out"``
#### NOTE: It's crucial to wrap the wildcard-containing path in quotes!

#### NOTE: This file will reject any file that contains the string "slurm."

#### Corin Wagen and Eugene Kwan, 2019

results = cctk.Ensemble()

parser = argparse.ArgumentParser(prog="analyze.py")
parser.add_argument("--g", action="store_true")
parser.add_argument("--h", action="store_true")
parser.add_argument("filename")
args = vars(parser.parse_args(sys.argv[1:]))

print("\n\033[3mreading files:\033[0m")

for filename in tqdm(sorted(glob.glob(args["filename"], recursive=True))):
    if re.search("slurm", filename):
        continue

    output_file = cctk.OrcaFile.read_file(filename)
    if output_file is None:
        continue
    if isinstance(output_file, list):
        output_file = output_file[-1]

    molecule = None
    if output_file.successful_terminations:
        results.add_molecule(*list(output_file.ensemble.items())[-1])
        molecule = output_file.ensemble.molecules[-1]
    elif len(output_file.ensemble) > 1:
        results.add_molecule(*list(output_file.ensemble.items())[-2])
        molecule = output_file.ensemble.molecules[-2]
    else:
        results.add_molecule(*list(output_file.ensemble.items())[-1])
        molecule = output_file.ensemble.molecules[-1]
    results[molecule, "iters"] = len(output_file.ensemble)
    results[molecule, "success"] = output_file.successful_terminations
    results[molecule, "imaginary"] = output_file.imaginaries()

if len(results) == 0:
    print("no jobs to analyze!")
    exit()

print("\n\n\033[3manalysis:\033[0m\n")
property_names = ["filename", "iters", "energy", "enthalpy", "quasiharmonic_gibbs_free_energy", "rms_step", "rms_gradient", "success", "imaginary"]
values = results[:, property_names]
if not isinstance(values[0], list):
    values = [values]

df = pd.DataFrame(values, columns=property_names).fillna("")

if args["g"]:
    df["rel_energy"] = (df.quasiharmonic_gibbs_free_energy- df.quasiharmonic_gibbs_free_energy.min()) * 627.509469
elif args["h"]:
    df["rel_energy"] = (df.enthalpy - df.enthalpy.min()) * 627.509469
else:
    df["rel_energy"] = (df.energy - df.energy.min()) * 627.509469

df.rename(columns={"rms_displacement": "rms_disp", "quasiharmonic_gibbs_free_energy": "GFE (corrected)"}, inplace=True)
df["filename"] = df["filename"].apply(lambda x: x[-60:])
df["GFE (corrected)"] = df["GFE (corrected)"].apply(lambda x: f"{x:.5f}")
df["rms_step"] = df["rms_step"].apply(lambda x: f"\033[92m{x}\033[0m" if float(x or 0) < 0.0001 else f"\033[93m{x}\033[0m")
df["rms_gradient"] = df["rms_gradient"].apply(lambda x: f"\033[92m{x}\033[0m" if float(x or 0) < 0.003 else f"\033[93m{x}\033[0m")
df["success"] = df["success"].apply(lambda x: f"\033[92m{x}\033[0m" if x else f"\033[93m{x}\033[0m")

df.columns = [f"\033[1m{c}\033[0m" for c in df.columns]
print(tabulate(df, headers="keys", tablefmt="presto", floatfmt=".5f"))

