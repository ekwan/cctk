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
        attach_to (int): atom number to replace with larger fragment
    """
    def __init__(self):
        pass    
