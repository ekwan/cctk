import sys
import os
import numpy as np

sys.path.append(os.path.relpath('../cctk'))

from cctk import XYZData, Molecule

output_file = XYZData.read_xyz('cctk/scripts/CpG.xyz')
output_file.title = "methylated CpG dinucleotide"

print(f"distance between atoms 1 and 2: {output_file.get_distance(1,2):.3f} A")

output_file.write_file('cctk/scripts/CpG_copy.xyz')
