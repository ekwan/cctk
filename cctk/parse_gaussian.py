import numpy as np
import re

from cctk.helper_functions import get_symbol

"""
Functions to help with parsing Gaussian files
"""


def read_geometries_and_energies(lines):
    """
    A large and unwieldy method - reads geometries, symbol lists, and energies from the file.

    Args:
        lines (list): the list of lines in the file

    Returns:
        array of geometries (each of which itself is an array of arrays)
        list of atomic symbols
        array of energies
        array containing the number of SCF iterations per step
    """

    file_geometries = []
    file_symbol_lists = []
    file_energies = []
    file_scf_iterations = []

    this_geometry = []
    this_symbol_list = []
    this_energy = None

    i = 0
    in_geometry_block = False
    while i < len(lines):
        # read the current line
        line = lines[i].strip()

        # detect geometry block
        if (line == "Standard orientation:") or (re.search("Cartesian Coordinates", line) or re.search("Input orientation:", line)):
            i += 5
            in_geometry_block = True
            continue
        elif in_geometry_block and line.startswith("------"):
            i += 1
            in_geometry_block = False
            continue

        # read geometry if applicable
        if in_geometry_block:
            fields = re.split(" +", line)
            if len(fields) == 6:
                if len(this_geometry) > 0 and fields[0] == "1":
                    # reset fields
                    this_geometry = []
                    this_symbol_list = []
                    this_energy = None
                try:
                    x, y, z = float(fields[3]), float(fields[4]), float(fields[5])
                    this_geometry.append([x, y, z])
                    symbol = get_symbol(fields[1])
                    this_symbol_list.append(symbol)
                except:
                    print("error parsing >>> %s" % line)
                    print(fields)
                    raise ValueError("error parsing geometry")

            elif len(fields) == 5:
                if len(this_geometry) > 0 and fields[0] == "1":
                    # reset fields
                    this_geometry = []
                    this_symbol_list = []
                    this_energy = None
                try:
                    x, y, z = float(fields[2]), float(fields[3]), float(fields[4])
                    this_geometry.append([x, y, z])
                    symbol = get_symbol(fields[1])
                    this_symbol_list.append(symbol)
                except:
                    print("error parsing >>> %s" % line)
                    print(fields)
                    raise ValueError("error parsing geometry")

            else:
                print("error parsing >>> " + line)
                raise ValueError("unexpected number of fields on geometry line")
        # read energy if applicable
        if not in_geometry_block and line.startswith("SCF Done"):
            fields = re.split(" +", line)
            if len(fields) != 9:
                print("error parsing >>> " + line)
                raise ValueError("unexpected number of fields on energy line")
            this_energy = float(fields[4])
            num_cycles = int(fields[7])

            if len(this_geometry) == 0:
                raise ValueError("energy without geometry")

            # store results
            file_geometries.append(this_geometry)
            file_symbol_lists.append(this_symbol_list)
            file_energies.append(this_energy)
            file_scf_iterations.append(num_cycles)

            # reinitialize arrays
            this_geometry = []
            this_symbol_list = []
            this_energy = None

        # go to next line
        i += 1

    # return result
    if len(file_symbol_lists) > 0:
        return file_geometries, file_symbol_lists, file_energies, file_scf_iterations
    else:
        (geometry, symbol_list) = extract_initial_geometry(lines)
        return [geometry], symbol_list, [], []


def search_for_block(lines, start, end, count=1, join=""):
    """
    Search through a file (lines) and locate a block starting with "start" (inclusive) and ending with "end" (exclusive).

    Args:
        lines (list): a list of the lines in the file
        start (str): a pattern that matches the start of the block (can contain special characters)
        end (str): a pattern that matches the end of the block (can contain special characters)
        count (int): how many matches to search for
        join (str): spacer between lines

    Returns:
        a single match (str) if count == 1 or a list of matches (str) if count > 1.
    """
    assert isinstance(count, int), "count needs to be an integer"
    assert isinstance(join, str), "join needs to be a string"

    current_match = ""
    match = [None] * count

    start_pattern = re.compile(start)
    end_pattern = re.compile(end)

    index = 0
    for line in lines:
        if current_match:
            if end_pattern.search(line):
                match[index] = current_match
                current_match = None
                index += 1

                if index == count:
                    break
            else:
                current_match = current_match + join + line.lstrip()
        else:
            if start_pattern.search(line):
                current_match = line.lstrip()

    if count == 1:
        return match[0]
    else:
        return match


def read_bonds(lines):
    """
    Reads bonding information from the output file.

    Args:
        lines (list): the list of lines in the file

    Returns:
        array of atoms that are bonded to each other
    """

    bond_array = []
    current_array = []

    i = 0
    in_bonding_section = False
    while i < len(lines):
        # read the current line
        line = lines[i].strip()

        if in_bonding_section == False:
            if re.search(r"Initial Parameters", line):
                i += 5
                in_bonding_section = True
                continue
            else:
                i += 1
                continue
        else:
            if re.search(r"! R", line):
                pieces = list(filter(None, line.split(" ")))
                atoms = pieces[2].replace("R", "").replace("(", "").replace(")", "").split(",")
                try:
                    current_array.append([int(atoms[0]), int(atoms[1])])
                except:
                    raise ValueError(f"error parsing line {i} - can't extract atoms!")
                i += 1
                continue
            elif (re.search(r"! A", line)) or (re.search("----", line)):
                in_bonding_section = False
                bond_array.append(current_array)
                current_array = []
                i += 1
                continue
            else:
                raise ValueError(f"can't parse line {i} for bonding!")

    return bond_array


def find_parameter(lines, parameter, expected_length, which_field):
    """
    Helper method to search through the output file and find key forces and displacements.

    Args:
        lines (list): list of lines in file
        parameter (string): test to search for
        expected_length (int): how many fields there should be
        which_field (int): which field the parameter is (zero-indexed)
    Returns:
        a list of all the extracted values
    """

    if (not isinstance(expected_length, int)) or (not isinstance(which_field, int)):
        raise TypeError("expected_length and which_field must be type int!")

    if which_field >= expected_length:
        raise ValueError("can't expect a field after the last field!")

    matches = []
    pattern = False

    try:
        pattern = re.compile(parameter)
    except:
        raise ValueError("pattern {pattern} cannot be compiled as a regex; try again!")

    if pattern:
        for line in lines:
            if pattern.search(line):
                fields = re.split(" +", line)
                fields = list(filter(None, fields))

                if len(fields) == expected_length:
                    try:
                        matches.append(float(fields[which_field]))
                    except:
                        pass
        return matches


def extract_initial_geometry(lines):
    """
    Helper method to search through the output file and find the initial geometry/symbol list, in cases where no SCF iterations finished.

    Args:
        lines (list): list of lines in file

    Returns:
        initial geometry (list)
        symbol list (list)
    """
    initial_geom_block = search_for_block(lines, "Symbolic Z-matrix:", "^ $", join="\n")
    atom_lines = initial_geom_block.split("\n")[2:]
    geometry = [None] * len(atom_lines)
    symbols = [None] * len(atom_lines)

    for idx, line in enumerate(atom_lines):
        frags = line.split()
        assert len(frags) == 4, f"too many fragments (more than 4) on line {line}"
        geometry[idx] = frags[1:]
        symbols[idx] = frags[0]
    return geometry, symbols
