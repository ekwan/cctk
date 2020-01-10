import sys, argparse, re, glob
import numpy as np

from cctk import GaussianFile, ConformationalEnsemble

#### Usage: ``python resubmit.py --type all --perturb "path/to/output/*.out"``
#### NOTE: It's crucial to wrap the wildcard-containing path in quotes!

#### NOTE: This file will reject any file that contains the string "slurm."

#### Corin Wagen and Eugene Kwan, 2019

parser = argparse.ArgumentParser(prog="resubmit.py")
parser.add_argument("filename")
args = vars(parser.parse_args(sys.argv[1:]))

assert args["filename"], "Can't resubmit files without a filename!"

ensembles = []

for filename in glob.iglob(args["filename"], recursive=True):
    if re.search("slurm", filename):
        continue

    output_file = GaussianFile.read_file(filename)
    ensembles.append(output_file.molecules)

new_ensemble = ConformationalEnsemble.join_ensembles(ensembles)

for mol in new_ensemble.molecules:
    cc_dist = mol.get_distance(1, 8)
    newfile = f"pop_{int(round(cc_dist*1000))}.gjf"
    GaussianFile.write_molecule_to_file(newfile, mol, "#p pop=hirshfeld m062x/6-31g(d)", None)
    print(f"generating {newfile}...")

