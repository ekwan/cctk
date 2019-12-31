#### Tutorial 2: Performing a Conformational Search on a Dipeptide

Predicting the conformation of short peptides in solution is a challenging and unsolved problem. 

##### Step 1: Generate Conformations to Search (`generate_conformers.py`)

Although in principle there are a practically infinite number of distinct structures that can be generated from 31 atoms, 
in practice most of low-energy conformational space can be sampled 
by selecting a few key rotatable bonds and letting Gaussian's `opt` keyword do the rest. 

For the purposes of this study, we chose to rotate four key bonds:
1. Rotating around F2Ala #1's alpha carbon with respect to the amide
1. Rotating around F2Ala #2's alpha carbon with respect to the amide
1. Rotating the difluoromethyl group in F2Ala #1
1. Rotating the difluoromethyl group in F2Ala #2

We also chose to sample *cis*/*trans* isomerism in each amide, as well as 

To get our starting `Molecule` object, we read from an `.xyz` file. Since `.xyz` files don't contain connectivity information, 
we have to generate the bonds automatically:

```
output_file = XYZFile.read_file('Ac-F2Ala-F2Ala-OMe.xyz')
output_file.molecule.assign_connectivity()
```

The actual heavy lifting is done by the following code, which creates copies of the reference structure with the selected
dihedral angles set to new values. 

```
current_idx = 0
        for structure in structures:
            for theta in thetas:
                new_structures[current_idx] = copy.deepcopy(structure.set_dihedral(*angles[idx], theta, check_result=False))
                current_idx += 1
```

Finally, the script writes the resultant structures to `.gjf` files

##### Step 2: Analyze Low-Level Results and Resubmit at Higher Level (`extract_unique.py`)


##### Step 3: Analyze High-Level Results (`analyze_final.py`)

<img src='dipeptide_conf_search/lowest_energy_conformer.png' width=600>
