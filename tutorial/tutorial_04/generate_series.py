import sys, os, argparse, glob, re, copy
import numpy as np

from cctk import GaussianFile, Molecule, Group, Ensemble
from cctk.load_groups import load_group

parser = argparse.ArgumentParser(prog="hammett_swap.py")
parser.add_argument("filename")
args = vars(parser.parse_args(sys.argv[1:]))
assert args["filename"], "Can't resubmit files without a filename!"

#### read in output file
output_file = GaussianFile.read_file(args["filename"])

#### read in genecp footer
footer = ""
with open('footer', 'r') as file:
    footer = file.read()

#### define groups and atoms
groups = ["NMe2", "OMe", "Me", "CF3", "CO2Me", "CN", "NO2", "F", "Cl", "SF5"]
p_atoms = [19, 30, 41]

ensemble = Ensemble()
headers = []
args = []

for group_name in groups:
    print(f"adding {group_name}")
    mol = copy.deepcopy(output_file.get_molecule().assign_connectivity())
    group = load_group(group_name)

    for atom in p_atoms:
        print(f"    adding to atom {atom}")
        mol = Group.add_group_to_molecule(mol, group, atom)

    ensemble.add_molecule(mol)
    headers.append(output_file.header)
    args.append({"footer": footer})

#### write everything to the same file with Link1
output_file.write_ensemble_to_file("hammett_NiCO3L.gjf", ensemble, headers, args)
