import sys
import re
import numpy as np
import networkx as nx
import math

from cctk import Molecule

class Group(Molecule):
    """
    Represents a functional group within a molecule - for *in silico* group substitutions. 

    Note that a Group instance does not need to be missing atoms. Rather, the atom given by `attach_to` will be replaced wholesale by another molecule, and the bond distances scaled automatically. 
    
    Attributes:
        attach_to (int): atom number to replace with larger fragment. must have only one bond! (e.g. H in F3C-H)
        adjacent (int): atom number that will be bonded to new molecule. (e.g. C in F3C-H)
    """
    def __init__(self, attach_to, **kwargs):
        super().__init__(**kwargs)
        self.add_attachment_point(attach_to)

    @classmethod
    def new_from_molecule(cls, molecule, attach_to):
        """
        Convenient method to convert ``molecule`` to ``group`` directly.
        """
        group = Group(attach_to, atoms=molecule.atoms, geometry=molecule.geometry, bonds=molecule.bonds.edges())
        return group

    def add_attachment_point(self, attach_to):
        """
        Adds ``attach_to``and ``adjacent`` attributes to the instance.
        
        Automatically centers atom ``adjacent`` on the origin, to simplify downstream mathematics.
        """
        if (len(super().get_adjacent_atoms(attach_to)) != 1): 
            raise ValueError(f"atom {attach_to} is making too many bonds!")
        
        self.attach_to = attach_to    

        adjacent = super().get_adjacent_atoms(attach_to)
        assert (len(adjacent) == 1), "can't substitute an atom with more than one adjacent atom!"
        self.adjacent = adjacent[0]
        
        adj_v = super().get_vector(self.adjacent)
        super().translate_molecule(-adj_v)
