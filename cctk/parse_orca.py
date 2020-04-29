import numpy as np
import re

from cctk.helper_functions import get_symbol, get_number
from cctk import OneIndexedArray, LazyLineObject

"""
Functions to help with parsing Orca files
"""
def read_geometries(lines, num_to_find):
    atomic_numbers = []
    geometries = []

    geom_blocks = lines.search_for_block("CARTESIAN COORDINATES \(ANGSTROEM\)", "CARTESIAN COORDINATES (A.U.)", join="\n", count=num_to_find)
    if num_to_find == 1:
        geom_blocks = [geom_blocks]

    for block in geom_blocks:
        rows = block.split("\n")
        numbers = []
        geometry = []

        for line in rows[2:]:
            if len(line.strip()) == 0:
                continue

            pieces = list(filter(None, line.split(" ")))

            if len(pieces) == 4:
                if re.match("[0-9]", pieces[0]):
                    numbers.append(int(pieces[0]))
                else:
                    numbers.append(int(get_number(pieces[0])))
                geometry.append([float(pieces[1]), float(pieces[2]), float(pieces[3])])

        atomic_numbers.append(OneIndexedArray(numbers, dtype=np.int8))
        geometries.append(OneIndexedArray(geometry))

        assert len(atomic_numbers) == len(geometries)
        for zs in atomic_numbers:
            assert np.array_equiv(zs, atomic_numbers[0])
        return atomic_numbers[0], geometries

def read_energies(lines):
    energies = lines.find_parameter("FINAL SINGLE POINT ENERGY", 5, 4)
    return energies

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



