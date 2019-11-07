import sys
import os
import numpy as np

sys.path.append(os.path.relpath('../cctk'))

from cctk import GaussianData, Molecule, GaussianJob

output_file = GaussianData.read_opt('cctk/scripts/acetaldehyde.out')

