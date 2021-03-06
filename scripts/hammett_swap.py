import sys, argparse, glob, re,

from cctk import GaussianFile
from cctk.load_groups import load_group
from cctk.Group import add_group_to_molecule

#### This is a script to automatically replace groups in molecules with predefined "Hammett"-type groups. 
#### Usage: ``python hammett_swap.py -a 9 13 -g trifluoromethyl -o bis_CF3_TS.gjf TS.out``

#### Corin Wagen and Eugene Kwan, 2019

parser = argparse.ArgumentParser(prog="hammett_swap.py")
parser.add_argument("--group", "-g", type=str)
parser.add_argument("--atoms", "-a", nargs='+')
parser.add_argument("--output", "-o", type=str)
parser.add_argument("filename")
args = vars(parser.parse_args(sys.argv[1:]))

assert args["filename"], "Can't resubmit files without a filename!"
assert args["group"], "Can't resubmit files without a filename!"
assert args["atoms"], "Can't resubmit files without a filename!"
for filename in glob.iglob(args["filename"], recursive=True):

    output_file = GaussianFile.read_file(filename)
    group = load_group(args['group'])
    mol = output_file.get_molecule().assign_connectivity()

    for atom in args['atoms']:
        mol = add_group_to_molecule(mol, group, atom)

    output_file.molecules[-1] = mol

    newfile = filename.rsplit('/',1)[-1]
    newfile = re.sub(r"out$", "gjf", newfile)
    if args["output"]:
        newfile = args["output"]

    output_file.write_file(newfile)
    print(f"{filename} > {newfile}")
