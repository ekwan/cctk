import os
import sys
import pickle

from cctk import MOL2File, Group
from cctk.group_substitution import add_group_to_molecule

filenames = [
    "MeH.mol2",
    "EtH.mol2",
    "iPrH.mol2",
    "tBuH.mol2",
    "OH2.mol2",
    "OMeH.mol2",
    "NHAcH.mol2",
    "NH3.mol2",
    "NMe2H.mol2",
    "CF3H.mol2",
    "HCN.mol2",
    "HNO2.mol2",
    "HCO2Me.mol2",
]

names = [
    "methyl",
    "ethyl",
    "isopropyl",
    "tert-butyl",
    "hydroxy",
    "methoxy",
    "acetamido",
    "amino",
    "dimethylamino",
    "trifluoromethyl",
    "cyano",
    "nitro",
    "carboxylmethyl",
]

groups_to_save = {}
for filename, name in zip(filenames, names):
    mol = MOL2File.read_file(f"groups/mol2/{filename}").molecules[0]
    mol.assign_connectivity()

    #### every molecule is set so you need to attach to atom 2
    new_group = Group.new_from_molecule(attach_to=2, molecule=mol)
    groups_to_save[name] = new_group

file = open('saved_groups.obj', 'wb')
pickle.dump(groups_to_save, file)
file.close()
