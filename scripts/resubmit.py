import sys, argparse, re, glob

from cctk import GaussianFile, OrcaFile

#### This is a script to resubmit failed Gaussian and ORCA files.
#### Parameters:
#### ``--type, -t``: which jobs to resubmit 
####     "failed": will resubmit jobs with no successes
####     "all": will resubmit all jobs (default)
#### ``--perturb, -p``: whether or not to apply a random geometric perturbation to each job
#### ``--output, -o``: name of output file (don't use with multiple output files!)

#### Usage: ``python resubmit.py --type all --perturb "path/to/output/*.out"``
#### NOTE: It's crucial to wrap the wildcard-containing path in quotes!

#### NOTE: This file will reject any file that contains the string "slurm."

#### Corin Wagen and Eugene Kwan, 2019
#### update to accept ORCA output, Joe Gair 2023

parser = argparse.ArgumentParser(prog="resubmit.py")
parser.add_argument("--type", "-t", type=str)
parser.add_argument("--perturb", "-p", action="store_true")
parser.add_argument("--output", "-o", type=str)
parser.add_argument("filename")
args = vars(parser.parse_args(sys.argv[1:]))

assert args["filename"], "Can't resubmit files without a filename!"

for filename in glob.iglob(args["filename"], recursive=True):
    if re.search("slurm", filename):
        continue

    try:
        if isinstance(GaussianFile.read_file(filename), GaussianFile):
            output_file = GaussianFile.read_file(filename)
            newfile_tail = "gjf"
        elif isinstance(OrcaFile.read_file(filename), OrcaFile):
            output_file = OrcaFile.read_file(filename)
            newfile_tail = "inp"
        else: print(f"error: {filename} is not recognized as an OrcaFile or a GaussianFile")
        if isinstance(output_file, list):
            output_file = output_file[-1]
        if args["perturb"]:
            output_file.get_molecule().perturb()
        success = output_file.successful_terminations

        if ((success == 0) and (args["type"] == "failed")) or (args["type"] == "all") or (args["type"] is None):
            newfile = filename.rsplit('/',1)[-1]
            newfile = re.sub(r"out$", newfile_tail, newfile)

            if args["output"]:
                newfile = args["output"]

            output_file.write_file(newfile)
            print(f"{filename} > {newfile}")

    except Exception as e:
        print(f"can't read file {filename}!\n{e}")
