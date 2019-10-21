import sys
import os
import numpy as np

sys.path.append(os.path.relpath('../cctk'))

from cctk import GaussianOptOutputFile

output_file = GaussianOptOutputFile('cctk/scripts/acetaldehyde.out')

energies = output_file.energies
scf_iter = output_file.scf_iterations
#rms_displacements = output_file.rms_displacements()
#rms_forces = output_file.rms_forces()

if output_file.successful: 
    print("Optimization converged!")

print("{0:2s} {1:15s} {2:15s} {3:4s}".format("#", "Energy (Hartree)", "Rel Energy (kcal)", "SCF Iterations"))

min_energy = np.min(energies)
for i, energy in enumerate(energies):
    rel_energy = (energy - min_energy) * 627.509 #### convert to kcal
    print("{0:2d} {1:15.5f} {2:15.5f} {3:4d}".format(i+1, energy, rel_energy, scf_iter[i]))
