import sys
import re
import glob
import numpy as np
import matplotlib.pyplot as plt

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

    idx1 = int((lg_dist - 1.5) * 10)
    idx2 = int((nu_dist - 1.2) * 10)

    print(f"{idx1} {idx2} {energy:.3f}")
    plot[idx1][idx2] = energy

min_energy = np.min(plot)
plot = plot - min_energy

fig = plt.figure(figsize=(10,7))
ax = fig.gca()

plt.imshow(plot, vmax=100, vmin=0)

#ax.set_ylim(bottom=1.5, top=3.3)
#ax.set_xlim(left=1.2, right=3.0)
#plt.savefig("plot.png")
plt.show()
