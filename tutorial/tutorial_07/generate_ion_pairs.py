import copy
import numpy as np
from cctk import Molecule, GaussianFile

num_structures = 25

footer = ""
with open('footer', 'r') as file:
    footer = file.read()

def spherical_random(radius=1):
    """
    Generates a random point on a sphere of radius ``radius``.
    """
    v = np.random.normal(size=3)
    v = v/np.linalg.norm(v)
    return v * radius

cation = GaussianFile.read_file("CuII-tBuBox-dication.out").get_molecule().center()
anion1 = GaussianFile.read_file("SbF6_anion.out").get_molecule().center()
anion2 = GaussianFile.read_file("OTf_anion.out").get_molecule().center()

anions = [anion1, anion2]
anion_names = ["SbF6", "OTf"]

for i in range(num_structures):
    trans_v = spherical_random(radius=8)
    for j in range(len(anions)):

        x = copy.deepcopy(anions[j])
        x.translate_molecule(trans_v)
        x.rotate_molecule(np.array([1,0,0]), np.random.random()*360)
        x.rotate_molecule(np.array([0,1,0]), np.random.random()*360)
        x.rotate_molecule(np.array([0,0,1]), np.random.random()*360)

        atoms = np.hstack((cation.atomic_numbers.T, x.atomic_numbers.T)).view(cctk.OneIndexedArray)
        geoms = np.vstack((cation.geometry, x.geometry)).view(cctk.OneIndexedArray)

        mx = Molecule(atomic_numbers=atoms, geometry=geoms, charge=1, multiplicity=2)
        GaussianFile.write_molecule_to_file(f"CuII-tBuBox-{anion_names[j]}_c{i}.gjf", mx, "#p opt b3lyp/genecp empiricaldispersion=gd3bj scrf=(smd, solvent=dichloromethane)", footer)
