import sys, os, argparse, glob, re, copy
import numpy as np

from cctk import GaussianFile, Molecule, Group, XYZFile
from cctk.group_substitution import add_group_to_molecule
from cctk.load_groups import load_group

parser = argparse.ArgumentParser(prog="hammett_swap.py")
parser.add_argument("--group", "-g", type=str)
parser.add_argument("--atoms", "-a", nargs='+')
parser.add_argument("filename")
args = vars(parser.parse_args(sys.argv[1:]))

assert args["filename"], "Can't resubmit files without a filename!"
assert args["group"], "Can't resubmit files without a filename!"
assert args["atoms"], "Can't resubmit files without a filename!"
for filename in glob.iglob(args["filename"], recursive=True):
 
    output_file = GaussianFile.read_file(filename)
    group = load_group(args['group'])

    print(output_file.get_molecule().atomic_numbers)
    mol = output_file.get_molecule()

    for atom in args['atoms']:
        print(atom)
        mol = add_group_to_molecule(mol, group, atom)
        print(mol.atomic_numbers)

    output_file.molecules[-1] = mol 
    print(output_file.get_molecule().atomic_numbers)

    newfile = filename.rsplit('/',1)[1]
    newfile = re.sub(r"out$", "gjf", newfile)
    output_file.write_file(newfile)
    print(f"{filename} > {newfile}")
