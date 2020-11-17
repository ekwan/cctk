try:
    import importlib.resources as pkg_resources
except ImportError:
    import importlib_resources as pkg_resources

from cctk import MOL2File, Group
from . import groups

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
    "FH.mol2",
    "ClH.mol2",
    "BrH.mol2",
    "IH.mol2",
    "SF5H.mol2",
    "SO3HH.mol2",
    "AcH.mol2",
    "CHOH.mol2",
]

names = [
    ["methyl", "Me", "CH3",],
    ["ethyl", "Et", "C2H5",],
    ["isopropyl", "iPr", "iC3H7",],
    ["tert-butyl", "tBu", "tC4H9",],
    ["hydroxy", "OH",],
    ["methoxy", "MeO", "OMe", "CH3O",],
    ["acetamido", "NHAc",],
    ["amino", "NH2",],
    ["dimethylamino", "Me2N", "NMe2",],
    ["trifluoromethyl", "CF3",],
    ["cyano", "CN",],
    ["nitro", "NO2",],
    ["carboxylmethyl", "MeO2C", "CO2Me",],
    ["fluoro", "F",],
    ["chloro", "Cl",],
    ["bromo", "Br",],
    ["iodo", "I",],
    ["pentafluorosulfanyl", "SF5",],
    ["sulfonyl", "SO3H",],
    ["acetyl", "Ac", "COMe",],
    ["formyl", "CHO",],
]

isomorphic = [
    [[3, 4, 5]],
    None,
    [[4, 8], [9, 10, 11, 5, 6, 7]],
    [[3, 7, 11], [4, 5, 6, 8, 9, 10, 12, 13, 14]],
    None,
    None,
    None,
    None,
    None,
    None,
    None,
    None,
    None,
    None,
    None,
    None,
    None,
    None,
    None,
    None,
    None,
]

def load_group(name):
    filename = None
    iso = None

    for row in names:
        if name in row:
            filename = filenames[names.index(row)]
            iso = isomorphic[names.index(row)]
            break

    assert filename is not None, f"can't find name {name}!"

    with pkg_resources.path(groups, filename) as file:
        mol = MOL2File.read_file(file).ensemble.molecules[0]
        mol.assign_connectivity()

        #### every molecule is set so you need to attach to atom 2
        new_group = Group.new_from_molecule(attach_to=2, molecule=mol, isomorphic=iso)
        return new_group

def group_iterator(symmetric_only=False):
    """
    Returns a generator over all *cctk*-predefined groups.
    """
    for row, iso in zip(names, isomorphic):
        if symmetric_only:
            if iso is None:
                continue
        yield load_group(row[0])
