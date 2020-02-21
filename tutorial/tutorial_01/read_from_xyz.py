import sys, argparse, re
from cctk import GaussianFile, XYZFile

#### Usage: python read_from_xyz.py -h "#p opt freq=noraman b3lyp/6-31g(d)" path/to/file.xyz

parser = argparse.ArgumentParser(prog="resubmit.py")
parser.add_argument("--header", "-h", type=str)
parser.add_argument("filename")
args = vars(parser.parse_args(sys.argv[1:]))

assert args["filename"], "Can't read file without a filename!"
assert args["header"], "Can't write file without a header!"

file = XYZFile.read_file(args["filename"])
newfile = args["filename"].rsplit('/',1)[-1]
newfile = re.sub(r"xyz$", "gjf", newfile)

GaussianFile.write_molecule_to_file(
    newfile,
    file.molecule,
    args["header"],
    None,
)
