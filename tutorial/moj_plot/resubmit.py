import sys, argparse
import re
import glob
import numpy as np

from cctk import GaussianFile, Molecule

#### This is a script to resubmit failed Gaussian files.
#### Parameters:
#### ``--type, -t``: which jobs to resubmit 
####     "failed": will resubmit jobs with no successes
####     "all": will resubmit all jobs
#### ``--perturb, -p``: whether or not to apply a random geometric perturbation to each job
#### ``--output, -o``: output file (only one)

#### Usage: ``python resubmit.py --type all --perturb "path/to/output/*.out"``
#### NOTE: It's crucial to wrap the wildcard-containing path in quotes!

#### NOTE: This file will reject any file that contains the string "slurm."

#### Corin Wagen and Eugene Kwan, 2019

parser = argparse.ArgumentParser(prog="resubmit.py")
parser.add_argument("--type", "-t", type=str)
parser.add_argument("--perturb", "-p", action="store_true")
parser.add_argument("--output", "-o", type=str)
parser.add_argument("filename")
args = vars(parser.parse_args(sys.argv[1:]))

cent = 1
lg = 7
nu = 8

assert args["filename"], "Can't resubmit files without a filename!"

for filename in glob.iglob(args["filename"], recursive=True):
    if re.search("slurm", filename):
        continue
    
    try: 
        output_file = GaussianFile.read_file(filename)
        if args["perturb"]:
            output_file.get_molecule().assign_connectivity()
            
            lg_dist = output_file.get_molecule().get_distance(cent, lg)
            nu_dist = output_file.get_molecule().get_distance(cent, nu)
            
            output_file.get_molecule().perturb()
            output_file.get_molecule().set_distance(cent, lg, lg_dist)
            output_file.get_molecule().set_distance(cent, nu, nu_dist)
        success = output_file.success

        if ((success == 0) and (args["type"] == "failed")) or (args["type"] == "all") or (args["type"] is None):
            newfile = filename.rsplit('/',1)[-1]
            newfile = re.sub(r"out$", "gjf", newfile)
            
            if args["output"]:
                newfile = args["output"]
            
            output_file.write_file(newfile)
            print(f"{filename} > {newfile}")
    except:
        print(f"Error parsing {filename} - will carry on regardless.")
