import sys, re, glob
import numpy as np

from cctk import GaussianFile, Molecule

#### This is a script to monitor the output of Gaussian files. 
#### In contrast to ``monitor.py``, this script analyzes many files! 
#### If the file has not successfully achieved SCF convergence at least once, that file will not display any information. 

#### Usage: ``python analyze.py "path/to/output/*.out"``
#### NOTE: It's crucial to wrap the wildcard-containing path in quotes!

#### NOTE: This file will reject any file that contains the string "slurm."

#### Corin Wagen and Eugene Kwan, 2019

filenames = sys.argv[1]
info = []
text_width = 60

for filename in sorted(glob.glob(filenames, recursive=True)):
    if re.search("slurm", filename):
        continue

    try:
        output_file = GaussianFile.read_file(filename)

        energy = output_file.energies[-1]
        gibbs = "          --  "
        try:
            gibbs = f"{output_file.gibbs_free_energy:16.4f}"
        except:
            pass
        iters = len(output_file.energies)
        rms_disp  = 0
        rms_force = 0

        try:
            rms_disp = output_file.rms_displacements[-1]
            rms_force = output_file.rms_forces[-1]
        except:
            pass

        success = "NO"
        if output_file.success:
            success = output_file.success

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

        info.append([filename[-text_width:], energy, energy * 627.509, iters, rms_force, rms_disp, success, imaginaries, gibbs])

    except:
        info.append([filename[-text_width:], 0, 0, 0, '', '', "NO", '', ''])


if len(info) > 0:
    min_energy = np.min([x[2] for x in info])
    def adjust_energy(row):
        if row[2] < 0:
            row[2] = row[2] - min_energy
        return row

    info = list(map(adjust_energy, info))

    print("{0:{text_width}}    {1:14}    {8:14}    {2:15}    {3:10}    {4:9}    {5:16}    {6:10}     {7:>12}".format(
        "File", "Energy (Hartree)", "Rel Energy (kcal)", "Iterations", "RMS Force", "RMS Displacement", "Success?", "Imaginaries?", "Gibbs (Hartree)", text_width=text_width
    ))

    for row in info:
        print("{0:{text_width}}    {1:14.4f}    {8:14s}    {2:15.2f}    {3:10}    {4:9}    {5:16}    {6:>10}     {7:>12}".format(*row, text_width=text_width))

else:
    print("no jobs to analyze!")
