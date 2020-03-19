import sys, re, glob
import numpy as np

from cctk import GaussianFile, Molecule
from cctk import parse_gaussian as parse

filenames = sys.argv[1]
info = []

for filename in sorted(glob.glob(filenames, recursive=True)):
    try:
        output_file, lines = GaussianFile.read_file(filename, return_lines=True)
    except:
        continue

    success = "NO"
    if output_file.success:
        success = output_file.success
    else:
        continue

    energy = output_file.energies[-1]
    iters = len(output_file.energies)
    mol = output_file.get_molecule()

    cation_anion_dist = 0
    mul_q = 0
    hir_q = 0
    if 29 in mol.atomic_numbers:
        mul_q = float(parse.find_parameter(lines, "    52  Cu", 4, 2)[-1])
        hir_q = float(parse.find_parameter(lines, "    52  Cu", 8, 2)[-1])

        if 16 in mol.atomic_numbers:
            cation_anion_dist = mol.get_distance(52, 54)
        elif 51 in mol.atomic_numbers:
            cation_anion_dist = mol.get_distance(52, 53)

    imaginaries = "--"
    try:
        if output_file.num_imaginaries() > 0:
            if output_file.num_imaginaries() > 1:
                imaginaries = ", ".join(output_file.imaginaries())
            else:
                imaginaries = output_file.imaginaries()[0]
    except:
        #### Will raise ValueError if job is not of type "FREQ"
        pass

    info.append([filename, energy, energy * 627.509, iters, mul_q, hir_q, success, imaginaries, cation_anion_dist])


if len(info) > 0:
    min_energy = np.min([x[2] for x in info])
    def adjust_energy(row):
        if row[2] < 0:
            row[2] = row[2] - min_energy
        return row

    info = list(map(adjust_energy, info))

    print("{0},{1},{2},{3},{4},{5},{6},{7},{8}".format(
        "File", "Energy (Hartree)", "Rel Energy (kcal)", "Iterations", "Mulliken Q", "Hirshfeld Q", "Success?", "Imaginaries?","X-Cu Distance"
    ))

    for row in info:
        print("{0},{1:.4f},{2:.2f},{3},{4:.4f},{5:.4f},{6},{7:},{8:.2f}".format(*row))

