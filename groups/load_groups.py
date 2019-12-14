import os
import sys
import pickle

try:
    import importlib.resources as pkg_resources
except ImportError:
    import importlib_resources as pkg_resources

from cctk import MOL2File, Group
from cctk.group_substitution import add_group_to_molecule

from ..groups import mol2  # relative-import the *package* containing the templates

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

def load_group(name):
    assert name in names, f"can't find group {name}!"
    filename = "groups/" + filnames[names.index(name)]

    isotope_file = pkg_resources.open_text(data, "isotopes.csv")
    mol = MOL2File.read_file(f"groups/mol2/{filename}").molecules[0]
    mol.assign_connectivity()

    #### every molecule is set so you need to attach to atom 2
    new_group = Group.new_from_molecule(attach_to=2, molecule=mol)
    return new_group

def pickle_groups():
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
