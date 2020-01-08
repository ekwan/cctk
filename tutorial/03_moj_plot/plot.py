import sys, re, glob
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter

from cctk import GaussianFile, Molecule, ConformationalEnsemble

#### Usage: ``python analyze.py "path/to/output/*.out"``
#### NOTE: It's crucial to wrap the wildcard-containing path in quotes!

#### NOTE: This file will reject any file that contains the string "slurm."

#### Corin Wagen and Eugene Kwan, 2019

filenames = sys.argv[1]

cent = 1
lg = 7
nu = 8

plot = np.zeros(shape=(18,18))

for filename in sorted(glob.glob(filenames, recursive=True)):
    if re.search("slurm", filename):
        continue

    output_file = GaussianFile.read_file(filename)
    energy = float(output_file.energies[-1]) * 627.509
    mol = output_file.get_molecule()

    lg_dist =  mol.get_distance(cent, lg)
    nu_dist =  mol.get_distance(cent, nu)

    idx1 = int(round((lg_dist - 1.5) * 10))
    idx2 = int(round((nu_dist - 1.2) * 10))
    plot[idx1][idx2] = energy

min_energy = np.min(plot)
plot = plot - min_energy

#### now to generate the graph...

fig = plt.figure(figsize=(10,8))
ax = fig.gca()

plt.imshow(plot, cmap="magma", vmax=60, vmin=0)

plt.title("More O’Ferrall–Jencks Plot for Nucleophilic Acyl Substitution", fontweight="bold")
ax.set_xlabel("Forming Bond (C1–C8)", fontweight="bold")
ax.set_ylabel("Breaking Bond (C1–Cl7)", fontweight="bold")

ax.set_xticks(range(18))
ax.set_yticks(range(18))
ax.set_xticklabels([f"{x:.2f}" for x in np.arange(1.2, 3.0, 0.1)])
ax.set_yticklabels([f"{y:.2f}" for y in np.arange(1.5, 3.3, 0.1)])

ax.set_xlim(left=0.5, right=17.5)
ax.set_ylim(top=0.5, bottom=17.5)

plt.colorbar(ax=ax).set_label("kcal/mol")

#plt.show()
plt.savefig("plot.png")
