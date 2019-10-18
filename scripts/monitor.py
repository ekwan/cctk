from cctk import GaussianOutputFile

output_file = GaussianOptOutputFile('file.out')

energies = output_file.energies()
rms_displacements = output_file.rms_displacements()
rms_forces = output_file.rms_forces()

if output_file.done(): 
    print("Optimization converged!")

for i, energy in enumerate(energies):
    print("{0:d} {1:.2f} {2:.2f} {3:.2f}".format(i, energy, rms_displacements[i], rms_forces[i])

