import sys, re, glob, cctk, argparse
import numpy as np
import pandas as pd
import multiprocessing as mp
from multiprocessing import freeze_support

from tabulate import tabulate
from tqdm import tqdm

#### This is a script to monitor the output of Gaussian files. 
#### In contrast to ``monitor.py``, this script analyzes many files! 
#### If the file has not successfully achieved SCF convergence at least once, that file will not display any information. 

#### Usage: ``python analyze.py "path/to/output/*.out"``
#### NOTE: It's crucial to wrap the wildcard-containing path in quotes!

#### NOTE: This file will reject any file that contains the string "slurm."

#### Corin Wagen and Eugene Kwan, 2019

def read(file):
    if re.match("slurm", file):
        return
    try:
        output_file = cctk.GaussianFile.read_file(file)
        if isinstance(output_file, list):
            if len(output_file) > 0:
                output_file = output_file[-1]
            else:
                return
        return output_file
    except Exception as e:
        print(f"Error reading {file}\n{e}")
        return

def main():
    results = cctk.Ensemble()

    parser = argparse.ArgumentParser(prog="analyze.py")
    parser.add_argument("--g", action="store_true")
    parser.add_argument("--h", action="store_true")
    parser.add_argument("--standard_state", action="store_true", default=False)
    parser.add_argument("--correct_gibbs", action="store_true", default=False)
    parser.add_argument("filename")
    args = vars(parser.parse_args(sys.argv[1:]))

    print("\n\033[3mreading files:\033[0m")

    files = glob.glob(args["filename"], recursive=True)

    pool = mp.Pool(processes=16)
    for output_file in tqdm(pool.imap(read, files), total=len(files)):
        molecule = None
        if output_file is None or (isinstance(output_file, list) and len(output_file) == 0):
            continue
        if output_file.successful_terminations:
            results.add_molecule(*list(output_file.ensemble.items())[-1])
            molecule = output_file.ensemble.molecules[-1]
        elif len(output_file.ensemble) > 1:
            results.add_molecule(*list(output_file.ensemble.items())[-2])
            molecule = output_file.ensemble.molecules[-2]
        elif len(output_file.ensemble) == 1:
            results.add_molecule(*list(output_file.ensemble.items())[-1])
            molecule = output_file.ensemble.molecules[-1]
        else:
            continue
        results[molecule, "iters"] = len(output_file.ensemble)
        results[molecule, "success"] = output_file.successful_terminations
        results[molecule, "imaginary"] = output_file.imaginaries()

    if len(results) == 0:
        print("no jobs to analyze!")
        exit()

    print("\n\n\033[3manalysis:\033[0m\n")
    property_names = None
    if args["correct_gibbs"]:
        property_names = ["filename", "iters", "energy", "enthalpy", "quasiharmonic_gibbs_free_energy", "rms_force", "rms_displacement", "success", "imaginary"]
    else:
        property_names = ["filename", "iters", "energy", "enthalpy", "gibbs_free_energy", "rms_force", "rms_displacement", "success", "imaginary"]

    values = results[:, property_names]
    if not isinstance(values[0], list):
        values = [values]

    df = pd.DataFrame(values, columns=property_names).fillna("")
    df.sort_values("filename", inplace=True)

    if args["g"]:
        df["rel_energy"] = (df.quasiharmonic_gibbs_free_energy - df.quasiharmonic_gibbs_free_energy.min()) * 627.509469
    elif args["h"]:
        df["rel_energy"] = (df.enthalpy - df.enthalpy.min()) * 627.509469
    else:
        try:
            df["rel_energy"] = (df.energy - df.energy.min()) * 627.509469
        except:
            df["rel_energy"] = 0

    if args["correct_gibbs"]:
        df.rename(columns={"rms_displacement": "rms_disp", "quasiharmonic_gibbs_free_energy": "GFE"}, inplace=True)
    else:
        df.rename(columns={"rms_displacement": "rms_disp", "gibbs_free_energy": "GFE"}, inplace=True)

    if args["standard_state"]:
        # convert to 1 M standard state = 6.354 e.u. * 298 K / 4184 J/kcal
        df["GFE"] = df["GFE"].apply(lambda x: f"{x - (0.4526 / 627.509):.5f}" if isinstance(x, float) else "")
    else:
        df["GFE"] = df["GFE"].apply(lambda x: f"{x:.5f}" if isinstance(x, float) else "")

    df["filename"] = df["filename"].apply(lambda x: x[-60:])
    df["rms_force"] = df["rms_force"].apply(lambda x: f"\033[92m{x}\033[0m" if float(x or 0) < 0.0001 else f"\033[93m{x}\033[0m")
    df["rms_disp"] = df["rms_disp"].apply(lambda x: f"\033[92m{x}\033[0m" if float(x or 0) < 0.003 else f"\033[93m{x}\033[0m")
    df["success"] = df["success"].apply(lambda x: f"\033[92m{x}\033[0m" if x else f"\033[93m{x}\033[0m")

    df.columns = [f"\033[1m{c}\033[0m" for c in df.columns]
    print(tabulate(df, headers="keys", tablefmt="presto", floatfmt=".5f"))

if __name__ == '__main__':
    freeze_support()
    main()

