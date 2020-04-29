import numpy as np
import re

from cctk.helper_functions import get_symbol
from cctk import OneIndexedArray, LazyLineObject

"""
Functions to help with parsing Orca files
"""
def read_geometries(lines):
    pass

def read_energies(lines):
    pass

def split_multiple_inputs(filename):
    """
    Splits ``filename`` into blocks by searching for _________.

    Args:
        filename (str): path to file

    Returns:
        list of list of ``LazyLineObject`` by input section
    """
    output_blocks = []

    start_block = 0
    with open(filename, "r") as lines:
        for idx, line in enumerate(lines):
            if re.search("Entering Link 1", line): # this will never be true for an Orca file -- this is just a stopgap
                output_blocks.append(LazyLineObject(file=filename, start=start_block, end=idx))
                start_block = idx
    output_blocks.append(LazyLineObject(file=filename, start=start_block, end=idx))

    return output_blocks



