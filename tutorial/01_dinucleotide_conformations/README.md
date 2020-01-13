#### Tutorial 1: Generating New Conformations

Generating new conformations for complex molecules can frequently be time-consuming, particularly when one desires to look at multiple rotatable bonds. 
For butane, there are only three minima about the C2–C3 bond, but looking at the equivalent minima for octane results in 3<sup>5</sup>=243 conformations— 
far too many to generate manually!
Scripting the creation of new conformers with *cctk* can thus be a powerful time-saving tool. 

For this example, we chose to study the *N*5-methylated CpG dinucleotide ((*N*5-methylation of cytosine is a common repressive epigenetic marker). 
New rotamers about the P–O bonds can easily be generated through the use of `molecule.set_dihedral()`:

```
angles = [0, 60, 120, 180, 240, 240.4, 360]
for x in angles:
    for y in angles:
        output_file.molecule.set_dihedral(1, 7, 6, 8, x)
        output_file.molecule.set_dihedral(23, 24, 25, 1, y)
        ensemble.add_molecule(output_file.molecule)
```

After elimination of redundant conformers, the resultant molecules are written to `.gjf` files:

```
for molecule in ensemble.molecules:
    molecule.check_for_conflicts()
    x = molecule.get_dihedral(1, 7, 6, 8)
    y = molecule.get_dihedral(23, 24, 25, 1)
    GaussianFile.write_molecule_to_file(f"conformers/CpG_{int(round(x))}_{int(round(y))}.gjf", molecule, "#p opt b3lyp/6-31g(d)", None)
```

The resultant files can then be submitted, and the resultant energies compared to determine the ground-state conformational distribution. 
