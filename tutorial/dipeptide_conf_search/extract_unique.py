import sys
import re
import glob
import numpy as np

from cctk import GaussianFile, Molecule, ConformationalEnsemble

#### This is a script to extract the lowest-energy unique conformers and resubmit them at a higher level of theory.
#### By default, this file will create new ``.gjf`` files for any conformer within 10 kcal/mol of the lowest-energy conformer. 

#### Usage: ``python extract_unique.py "path/to/output/*.out"``
#### NOTE: It's crucial to wrap the wildcard-containing path in quotes!

#### Corin Wagen and Eugene Kwan, 2019

filenames = sys.argv[1]
info = []
text_width = 70

to_rotate = [[1, 3, 5, 7], [9, 11, 13, 15], [5, 3, 6, 8], [12, 11, 14, 16]]

ensemble = ConformationalEnsemble()

for filename in sorted(glob.glob(filenames, recursive=True)):
    if re.search("slurm", filename):
        continue
   
    try:  
        output_file = GaussianFile.read_file(filename)
        
        if len(output_file.energies) > 0:
            mol = output_file.get_molecule() 
            ensemble.add_molecule(mol, energy=output_file.energies[-1]*627.509)
    except:
        print(f"skipping f{filename} due to error...")

print(f"{len(ensemble.molecules)} conformers before elimination of redundant")
ensemble.eliminate_redundant(cutoff=0.6)
print(f"{len(ensemble.molecules)} conformers after elimination of redundant")

best_confs = ensemble.get_within_cutoff(cutoff=10)
for idx, molecule in enumerate(best_confs):
    GaussianFile.write_molecule_to_file(f"conformer_v2_{idx:03d}.gjf", molecule, "#p opt freq=noraman m062x/6-31g(d) scrf=(smd,solvent=diethylether)", None)
    
for idx, molecule in enumerate(list(ensemble.molecules[np.argsort(ensemble.energies)])):
    items = [idx+1, np.sort(ensemble.energies)[idx]]    

    for atoms in to_rotate:
        items.append(molecule.get_dihedral(*atoms))

    if idx==0:
        print("Molecule    Energy    D" + str(" D".join(str(a) for a in to_rotate)))

    print(str("        ".join(str(x)[0:8] for x in items)))
