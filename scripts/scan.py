import cctk, yaml, tqdm, argparse, sys, itertools, os, tqdm

parser = argparse.ArgumentParser(prog="scan.py")
parser.add_argument("scanfile")
parser.add_argument("filenames", nargs="+")

args = vars(parser.parse_args(sys.argv[1:]))

settings = dict()
with open(args["scanfile"], "r+") as f:
    settings = yaml.safe_load(f)

assert "base" in settings
assert "dimensions" in settings

def build_file_with_coords(file, coords, name):
    route_card, footer = settings["base"]["route_card"], settings["base"]["footer"]

    name = os.path.splitext(name)[0]

    assert len(settings["dimensions"]) == len(coords)
    for dimension, coord in zip(settings["dimensions"].values(), coords):
         name += "_" + dimension["name"][coord]
         if "route_card" in dimension and dimension["route_card"][coord] is not None:
            if route_card:
                route_card += " " + dimension["route_card"][coord]
            else:
                route_card = dimension["route_card"][coord]
         if "footer" in dimension and dimension["footer"][coord]is not None:
            if footer:
                footer += " " + dimension["footer"][coord]
            else:
                footer = dimension["footer"][coord]

    file.write_file(name + ".gjf", route_card=route_card, footer=footer)

for filename in args["filenames"]:
    file = cctk.GaussianFile.read_file(filename)
    for coords in tqdm.tqdm(itertools.product(*[range(len(dim["name"])) for dim in settings["dimensions"].values()])):
        build_file_with_coords(file, coords, filename)
