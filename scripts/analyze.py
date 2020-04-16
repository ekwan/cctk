import sys, re, glob, cctk
import numpy as np

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
    for p in ["energy", "enthalpy", "gibbs_free_energy", "rms_force", "rms_displacement"]:
        x = output_file.ensemble[output_file.get_molecule(), p]
        if x is not None:
            try:
                results[output_file.get_molecule(), p] = f"{float(results[output_file.get_molecule(), p]):.5f}"
            except:
                results[output_file.get_molecule(), p] = results[output_file.get_molecule(), p]
        else:
            results[output_file.get_molecule(), p] = ""

print(f"\033[1m{'filename':50} {'iters':>5} {'E':>15} {'∆E':>11} {'G':>15} {'H':>15} {'force':>8} {'disp':>8} {'done?':>8} {'# imag':>8}\033[0m")
min_energy = min([float(x) for x in results[:,"energy"]])
for mol, row in results:
    rel_energy = f"{(float(row['energy']) - min_energy) * 627.509:8.2f}"
    print(f"{row['filename']:50} {row['iterations']:>5} {row['energy']:>15} {rel_energy:>11} {row['gibbs_free_energy']:>15} {row['enthalpy']:>15} {row['rms_force']:>8} {row['rms_displacement']:>8} {row['success']:>8} {row['num_imag']:>8}")

if len(results) == 0:
    print("no jobs to analyze!")
