import numpy as np
import re

from cctk.helper_functions import get_number
from cctk import OneIndexedArray, LazyLineObject

"""
Functions to help with parsing Orca files
"""
def read_geometries(lines, num_to_find):
    atomic_numbers = []
    geometries = []

    geom_blocks = lines.search_for_block("CARTESIAN COORDINATES \(ANGSTROEM\)", "CARTESIAN COORDINATES", join="\n", count=num_to_find, max_len=1000)
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
    iters = lines.find_parameter("SCF CONVERGED AFTER", 7, 4)
    return energies, iters

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

def read_mulliken_charges(lines):
    """
    Reads charges.

    Args:
        lines (list): list of lines in file

    Returns:
        ``cctk.OneIndexedArray`` of charges
    """
    charges = []
    charge_block = lines.search_for_block("MULLIKEN ATOMIC CHARGES", "Sum of atomic charges", join="\n")
    for line in charge_block.split("\n")[2:]:
        fields = re.split(" +", line)
        fields = list(filter(None, fields))

        if len(fields) == 4:
            charges.append(float(fields[3]))

    return OneIndexedArray(charges)


def read_loewdin_charges(lines):
    """
    Reads charges.

    Args:
        lines (list): list of lines in file

    Returns:
        ``cctk.OneIndexedArray`` of charges
    """
    charges = []
    charge_block = lines.search_for_block("LOEWDIN ATOMIC CHARGES", "^$", join="\n")
    for line in charge_block.split("\n")[2:]:
        fields = re.split(" +", line)
        fields = list(filter(None, fields))

        if len(fields) == 4:
            charges.append(float(fields[3]))

    return OneIndexedArray(charges)

def read_header(lines):
    for line in lines:
        if re.match("!", line):
            return line

def read_blocks_and_variables(lines):
    blocks = {}
    variables = {}

    current_key = None
    current_val = []
    for line in lines:
        if current_key is not None:
            if re.match("end", line):
                blocks[current_key] = current_val
                current_key = None
                current_val = []
            else:
                current_val.append(line)
                continue
        if re.match("%", line):
            fields = re.split(" +", line.lstrip("%"))
            if len(fields) == 1:
                current_key = fields[0]
            else:
                variables[fields[0]] = " ".join(fields[1:])

    return variables, blocks

def extract_input_file(lines):
    input_block = lines.search_for_block("INPUT FILE", "\*\*\*\*END OF INPUT\*\*\*\*", join="\n")
    input_lines = []
    for line in input_block.split("\n")[3:]:
        [_, line] = line.split(">")
        line = line.lstrip()
        input_lines.append(line)
    return input_lines

def read_freqs(lines):
    freq_block = lines.search_for_block("VIBRATIONAL FREQUENCIES", "NORMAL MODES", join="\n", max_len=1000)
    if freq_block is None:
        return []
    freqs = []
    for line in freq_block.split("\n"):
        fields = re.split(" +", line.strip())
        if len(fields) == 3:
            if fields[2] == "cm**-1" and float(fields[1]) > 0:
                freqs.append(float(fields[1]))
    return freqs

def read_gradients(lines, num_to_find):
    grad_blocks = lines.search_for_block("Geometry convergence", "Max\(Bonds", join="\n", count=num_to_find)
    if grad_blocks is None:
        return

    rms_grad = []
    max_grad = []
    rms_step = []
    max_step = []
    for grad_block in grad_blocks:
        if grad_block is None:
            continue
        for line in grad_block.split("\n"):
            fields = re.split(" +", line.strip())
            if len(fields) == 5:
                if fields[0] == "RMS" and fields[1] == "gradient":
                    rms_grad.append(float(fields[2]))
                if fields[0] == "MAX" and fields[1] == "gradient":
                    max_grad.append(float(fields[2]))
                if fields[0] == "RMS" and fields[1] == "step":
                    rms_step.append(float(fields[2]))
                if fields[0] == "MAX" and fields[1] == "step":
                    max_step.append(float(fields[2]))

    return rms_grad, max_grad, rms_step, max_step

def read_nmr_shifts(lines, num_atoms):
    """
    Helper method to search through output file and read NMR shifts.

    Args:
        lines (list): list of lines in file
        num_atoms (int): number of atoms expected

    Returns:
        list of isotropic NMR shifts (np.ndarray)
    """
    # assumes that lines only come from one Link1 section
    shieldings = []
    block = lines.search_for_block("Nucleus  Element", "^$", join="\n")
    for line in block.split("\n")[2:]:
        fields = line.split()
        if len(fields) == 4:
            try:
                shielding = float(fields[2])
            except:
                raise ValueError(f"Error parsing NMR shielding output line:\n{line}")
            shieldings.append(shielding)

    if len(shieldings) != 0:
        assert len(shieldings) == num_atoms, f"Expected {num_atoms} shieldings but found {len(shieldings)}!"
        return np.asarray(shieldings).view(OneIndexedArray)
    else:
        #### we can catch this problem later if the file is finished
        return None


