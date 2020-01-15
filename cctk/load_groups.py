import os
import sys

try:
    import importlib.resources as pkg_resources
except ImportError:
    import importlib_resources as pkg_resources

from cctk import MOL2File, Group
from cctk.group_substitution import add_group_to_molecule

from . import data, groups

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
    filename = filenames[names.index(name)]

    with pkg_resources.path(groups, filename) as file:
        mol = MOL2File.read_file(file).molecules[0]
        mol.assign_connectivity()

        #### every molecule is set so you need to attach to atom 2
        new_group = Group.new_from_molecule(attach_to=2, molecule=mol)
        return new_group
