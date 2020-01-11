import sys, re, glob
import numpy as np
import matplotlib.pyplot as plt

from cctk import GaussianFile, Molecule
import cctk.parse_gaussian as parse

#### Usage: ``python analyze_dipole.py "path/to/output/*.out"``
#### NOTE: It's crucial to wrap the wildcard-containing path in quotes!

#### NOTE: This file will reject any file that contains the string "slurm."

#### Corin Wagen and Eugene Kwan, 2019

filenames = sys.argv[1]

energies = {}
dipole = {}

for filename in sorted(glob.glob(filenames, recursive=True)):
    if re.search("slurm", filename):
        continue

    (output_file, lines) = GaussianFile.read_file(filename, return_lines=True)
    dist = int(round(output_file.get_molecule().get_distance(1, 8) * 1000))

    energies[dist] = output_file.energies[-1]

    #### look for the dipole moment
    try:
        dipole_line = parse.search_for_block(lines, "Dipole", "Quadrupole")

        fields = re.split(" +", dipole_line)
        fields = list(filter(None, fields))
        dipole[dist] = float(fields[-1])

    except:
        pass

#### generate graphs
min_energy = np.min(list(energies.values()))
energies = {k: (e - min_energy) * 627.509 for k, e in energies.items()}

fig = plt.figure(figsize=(10,7))
ax1 = fig.gca()
ax1.scatter(list(energies.keys()), list(energies.values()), c='black', alpha=0.8, label="Energy")
ax1.set_ylim(top=30, bottom=0)
ax1.set_xlabel("C1-C5 Distance (mÅ)")
ax1.set_ylabel("Energy (kcal/mol; M06-2X)")

ax2 = ax1.twinx()
ax2.scatter(list(dipole.keys()), list(dipole.values()), c='blue', alpha=0.8, label="Dipole")
ax2.set_ylim(top=3, bottom=0)
ax2.set_ylabel("Dipole Moment (M06-2X)")

plt.title("Change in Dipole Moment over IRC")
plt.legend(loc="best")

plt.show()
plt.savefig('dipole_graph.png')
