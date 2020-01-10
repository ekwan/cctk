import sys, re, glob
import numpy as np
import matplotlib.pyplot as plt

from cctk import GaussianFile, Molecule
import cctk.parse_gaussian as parse


#### Usage: ``python analyze.py "path/to/output/*.out"``
#### NOTE: It's crucial to wrap the wildcard-containing path in quotes!

#### NOTE: This file will reject any file that contains the string "slurm."

#### Corin Wagen and Eugene Kwan, 2019

filenames = sys.argv[1]
info = []
text_width = 70

energies = {}
nics = {}
C1_charge = {}
C5_charge = {}
O7_charge = {}
C8_charge = {}
C9_charge = {}

for filename in sorted(glob.glob(filenames, recursive=True)):
    if re.search("slurm", filename):
        continue

    (output_file, lines) = GaussianFile.read_file(filename, return_lines=True)
    dist = int(round(output_file.get_molecule().get_distance(1, 5) * 1000))

    energies[dist] = output_file.energies[-1]

    try:
        nics[dist] = parse.find_parameter(lines, "17  Bq   Isotropic", 8, 4)[0]
    except:
        pass

    try:
        C1_charge[dist] = parse.find_parameter(lines, "     1  C", 8, 2)[-1]
        C5_charge[dist] = parse.find_parameter(lines, "     5  C", 8, 2)[-1]
        O7_charge[dist] = parse.find_parameter(lines, "     7  O", 8, 2)[-1]
        C8_charge[dist] = parse.find_parameter(lines, "     8  C", 8, 2)[-1]
        C9_charge[dist] = parse.find_parameter(lines, "     9  C", 8, 2)[-1]
    except: 
        pass
     

