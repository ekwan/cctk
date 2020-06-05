import numpy as np
import re
import ahocorasick
import itertools
import string

import cctk
from cctk.helper_functions import get_symbol, get_number, get_corrected_free_energy

"""
Functions to help with parsing Gaussian files
"""

def read_geometries_and_energies(lines_obj):
    """
    A large and unwieldy method - reads geometries, symbol lists, and energies from the file.

    Args:
        lines (list): the list of lines in the file

    Returns:
        array of geometries (each of which itself is an array of arrays)
        list of atomic symbols
        array of energies
        array containing the number of SCF iterations per step
        num successful terminations
        elapsed time
    """

    # abandon all hope, ye who enter here

    file_geometries = []
    file_symbol_lists = []
    file_energies = []
    file_scf_iterations = []

    this_geometry = []
    this_symbol_list = []
    this_energy = None

    elapsed_time = 0
    success = 0
    i = 0

    #### we'll read this into memory to speed this up
    lines = list(lines_obj)

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
#        elif not in_geometry_block:
#            if line.startswith("Elapsed time"):
#                fields = line.split()
#                assert len(fields) == 10, f"unexpected number of fields on elapsed time line:\n{line}"
#                days = float(fields[2])
#                hours = float(fields[4])
#                minutes = float(fields[6])
#                seconds = float(fields[8])
#                elapsed_time += days * 86400 + hours * 3600 + minutes * 60 + seconds
#
#            elif line.startswith("Normal termination"):
#                success += 1

        # go to next line
        i += 1

    # return result
    if len(file_symbol_lists) > 0:
        return file_geometries, file_symbol_lists[0], file_energies, file_scf_iterations, success, elapsed_time
    else:
        (geometry, symbol_list) = extract_initial_geometry(lines_obj)
        return [geometry], symbol_list, [], [], success, elapsed_time

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

    #### we'll read this into memory to speed this up
    lines = list(lines)

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

    if len(bond_array):
        return bond_array[0]
    else:
        return None

def extract_initial_geometry(lines):
    """
    Helper method to search through the output file and find the initial geometry/symbol list, in cases where no SCF iterations finished.

    Args:
        lines (list): list of lines in file

    Returns:
        initial geometry (list)
        symbol list (list)
    """
    initial_geom_block = lines.search_for_block("Symbolic Z-matrix:", "^ $", join="\n", max_len=1000)

    if initial_geom_block is not None:
        atom_lines = initial_geom_block.split("\n")[2:]
        geometry = [None] * len(atom_lines)
        symbols = [None] * len(atom_lines)

        for idx, line in enumerate(atom_lines):
            frags = line.split()
            assert len(frags) == 4, f"too many fragments (more than 4) on line {line}"
            geometry[idx] = frags[1:]
            symbols[idx] = frags[0]
        return geometry, symbols

    else:
        initial_geom_block = lines.search_for_block("Standard orientation:", "Rotational", join="\n", max_len=1000)

        if initial_geom_block is None:
            raise ValueError("can't find valid geometry")

        atom_lines = initial_geom_block.split("\n")[5:-1] # gaussian, y u be like dis??
        geometry = [None] * len(atom_lines)
        symbols = [None] * len(atom_lines)

        for idx, line in enumerate(atom_lines):
            frags = line.split()
            assert len(frags) == 6, f"too many fragments (more than 6) on line {line}"
            geometry[idx] = frags[3:]
            symbols[idx] = frags[1]

            return geometry, [get_symbol(z) for z in symbols]

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
    for line in lines:
        fields = line.split()
        if len(fields) == 8 and fields[2] == "Isotropic" and fields[5] == "Anisotropy":
            fields = line.split()
            assert len(fields) == 8, f"Expected 8 fields on an NMR shielding output line but found {len(fields)} instead!"
            try:
                shielding = float(fields[4])
            except:
                raise ValueError(f"Error parsing NMR shielding output line:\n{line}")
            shieldings.append(shielding)

    if len(shieldings) is not 0:
        assert len(shieldings) == num_atoms, f"Expected {num_atoms} shieldings but found {len(shieldings)}!"
        return np.asarray(shieldings)
    else:
        #### we can catch this problem later if the file is finished
        return None

def split_link1(filename):
    """
    Splits ``filename`` into blocks by searching for "Entering Link 1".

    Args:
        filename (str): path to file

    Returns:
        list of list of lines by Link1 section; so a file with one Link1 specification would return [lines1, lines2]
    """
    link1_blocks = []

    start_block = 0
    with open(filename, "r") as lines:
        for idx, line in enumerate(lines):
            if re.search("Entering Link 1", line):
                link1_blocks.append(cctk.LazyLineObject(file=filename, start=start_block, end=idx))
                start_block = idx
    link1_blocks.append(cctk.LazyLineObject(file=filename, start=start_block, end=idx))

    return link1_blocks[1:] #### the first block is just a few lines

def extract_link0(lines):
    """
    Extracts Link 0 commands from ``lines``.

    Args:
        lines (list): list of lines in file

    Returns:
        dictionary of link 0 commands
    """
    output = {}

    #### read first 200 lines of file to memory
    lines = [x for _, x in zip(range(200), lines)]

    #### find where the header starts
    for idx, line in enumerate(lines):
        if re.match(" #p", line):
            break
    end_idx = idx - 1


    ##### then go back until the big header that says "GAUSSIAN 16"
    while not re.match(" \*\*\*\*\*\*\*\*", lines[idx]):
        idx += -1
        if idx < 10:
            raise ValueError("missed the termination point somehow")
    start_idx = idx

    #### and search within that little block for lines matching the right description
    for line in lines[start_idx:end_idx]:
        if re.match(" \%", line):
            pieces = line[2:].split("=")
            output[pieces[0]] = pieces[1]

    return output

def read_forces(lines):
    """
    Reads forces from a Gaussian ``force`` job.

    Args:
        lines (list): list of lines in file

    Returns:
        (n x 3) ``cctk.OneIndexedArray`` of forces
    """
    forces = []
    force_block = lines.search_for_block("Forces \(Hartrees/Bohr\)", "Cartesian Forces", join="\n", max_len=1000)
    for line in force_block.split("\n")[2:]:
        fields = re.split(" +", line)
        fields = list(filter(None, fields))

        if len(fields) == 5:
            forces.append([float(fields[2]), float(fields[3]), float(fields[4])])

    return cctk.OneIndexedArray(forces)

def read_mulliken_charges(lines):
    """
    Reads charges from a Gaussian ``pop`` job.

    Args:
        lines (list): list of lines in file

    Returns:
        ``cctk.OneIndexedArray`` of charges
    """
    charges = []
    charge_block = lines.search_for_block(" Mulliken charges:", " Sum of Mulliken charges", join="\n", max_len=1000)
    if charge_block is not None:
        for line in charge_block.split("\n")[2:]:
            fields = re.split(" +", line)
            fields = list(filter(None, fields))

            if len(fields) == 3:
                charges.append(float(fields[2]))

    return cctk.OneIndexedArray(charges)

def read_hirshfeld_charges(lines):
    """
    Reads charges from a Gaussian ``pop`` job.

    Args:
        lines (list): list of lines in file

    Returns:
        ``cctk.OneIndexedArray`` of charges
        ``cctk.OneIndexedArray`` of spin densities
    """
    charges = []
    spins = []
    charge_block = lines.search_for_block("Hirshfeld charges, spin densities, dipoles, and CM5 charges", " Hirshfeld charges", join="\n", max_len=1000)
    for line in charge_block.split("\n")[2:]:
        fields = re.split(" +", line)
        fields = list(filter(None, fields))

        if len(fields) == 8:
            charges.append(float(fields[2]))
            spins.append(float(fields[3]))

    return cctk.OneIndexedArray(charges), cctk.OneIndexedArray(spins)

def read_dipole_moment(lines):
    """
    Reads dipole moment from a Gaussian job.

    Args:
        lines (list): list of lines in file

    Returns:
        dipole moment (magnitude only)
    """
    dipole_block = lines.search_for_block(" Electronic spatial extent", " Quadrupole moment", join="\n")
    for line in dipole_block.split("\n")[1:]:
        fields = re.split(" +", line)
        fields = list(filter(None, fields))

        if len(fields) == 8:
            return float(fields[7])

def read_j_couplings(lines, n_atoms):
    """
    Helper method to search through output file and read J couplings

    Args:
        lines (list): list of lines in file
        n_atoms (int): how many atoms are in the molecule

    Returns:
        ``couplings`` symmetric 2D np.array of couplings (in Hz) with zero-indexed atoms on both axes
        or None if no couplings were found
    """
    couplings = np.zeros((n_atoms,n_atoms))
    n_full_blocks, lines_in_partial_block = divmod(n_atoms,5)
    n_lines = 5 * (n_full_blocks * (n_full_blocks+1) / 2) + n_full_blocks + 1
    if lines_in_partial_block > 0:
        n_lines += 1 + lines_in_partial_block
    n_lines = int(n_lines)

    #### Need to refactor this

    if isinstance(lines, cctk.LazyLineObject):
        lines = lines.search_for_block("Total nuclear spin-spin coupling J \(Hz\):", None, max_len=n_lines, join="\n")
        if lines is None:
            return None
        lines = lines.split("\n")

        assert n_lines == len(lines), f"found {len(lines)} lines but expected {n_lines} lines!"
    elif isinstance(lines, list):
        lines = lines[0].split("\n")

    i = 0
    read_column_indices = False
    read_row = False
    this_column_indices = []
    while i < n_lines:
        # get current line
        line = lines[i]

        # if this is the header, we should be reading the column indices next
        if "Total nuclear spin-spin coupling J (Hz):" in line:
            i += 1
            read_column_indices = True
            continue

        # this is not the header, so split the fields
        fields = line.split()

        # read the column indices
        if read_column_indices:
#            this_n_columns = len(fields)
            this_column_indices = [ int(j)-1 for j in fields ]
            i += 1
            read_column_indices = False
            read_row = True
            continue
        elif read_row:
            row = int(fields[0])-1
            for j,value in enumerate(fields[1:]):
                column = this_column_indices[j]
                value = value.replace("D","E")
                value = float(value)
                couplings[row,column] = value
                couplings[column,row] = value

            # check if we have read the entire matrix
            if row == n_atoms - 1 and column == n_atoms - 1:
                break

            # check if this is the end of the current block
            if row == n_atoms - 1:
                read_column_indices = True
                read_row = False
                i += 1
                continue

            read_row = True
            i += 1
            continue
        else:
            raise ValueError("impossible")

    return couplings

def extract_success_and_time(lines):
    elapsed_time = 0
    success = lines.find_parameter("Normal termination", 11, 0)
    time = lines.find_parameter("Elapsed time", 10, [2, 4, 6, 8])
    for fields in time:
        elapsed_time += fields[0] * 86400 + fields[1] * 3600 + fields[2] * 60 + fields[3]
    return len(success), elapsed_time

def read_file_fast(file_text, filename, link1idx, max_len=1000, extended_opt_info=False):

    #### "Make your bottleneck routines fast, everything else clear" - M. Scott Shell, UCSB
    #### Welcome to the fast part!

    #### Here we identify all the lines we're going to scrape
    words = [
        "SCF Done",
        "Entering Link 1",
        "Normal termination",
        "Elapsed time",
        "Multiplicity",
        "RMS     Force", #5
        "RMS     Displacement",
        "Maximum Force",
        "Maximum Displacement",
        "Cartesian Forces",
        "Internal  Forces", #10
        "Predicted change in Energy",
        "thermal Enthalpies",
        "thermal Free Energies",
        "Frequencies",
        "Temperature", #15
        "Isotropic",
    ]

    #### And here are the blocks of text
    #### format: [start, stop, num]

    blocks = [
        ["#p", "----", 1],
        ["l101.exe", "Symbolic Z-matrix", 1],
        ["The following ModRedundant input section", "\n \n", 1],
        [
            ["Input orientation", "Standard orientation", "Cartesian Coordinates"],
           "Rotational",
            1000,
        ],
        ["Wallingford", "#p", 1],
        ["Initial Parameters", "! A", 1], #5
        ["Total nuclear spin-spin coupling J", "Leave Link", 1],
        ["Forces (Hartrees/Bohr)", "Cartesian Forces", 1],
        ["Hirshfeld charges, spin densities, dipoles, and CM5 charges", " Hirshfeld charges", 1],
        ["Mulliken charges", "Sum of Mulliken charges", 1],
        ["Electronic spatial extent", "Quadrupole moment", 1], #10
    ]

    word_matches = [[] for _ in words]
    block_matches = [[] for _ in blocks]

    A = ahocorasick.Automaton()

    for idx, word in enumerate(words):
        A.add_word(word, idx)

    for idx, b in enumerate(blocks):
        if isinstance(b[0], list):
            for start in b[0]:
                A.add_word(start, ("start", idx))
        else:
            A.add_word(b[0], ("start", idx))

    #### perform search
    A.make_automaton()
    found_words = A.iter(file_text)

    #### now, we have to expand our one-character matches to whole lines/blocks
    #### this is the slowest part
    for position, idx in found_words:
        if isinstance(idx, int):
            stepsize = 10

            match = file_text[position]
            i = position + 1
            while match[-1-stepsize:].find("\n") < 0:
                match = match + file_text[i:i+stepsize]
                i += stepsize

            match = match.split("\n")[0]

            j = position
            while match[:stepsize].find("\n") < 0:
                match = file_text[j-stepsize:j] + match
                j += -1 * stepsize

            match = match.split("\n")[-1]
            word_matches[idx].append(match)

        elif isinstance(idx, tuple):
            idx = idx[1]
            if len(block_matches[idx]) >= blocks[idx][2]:
                continue

            match = ""
            i = position - len(blocks[idx][0]) + 1
            end = blocks[idx][1]

            stepsize = 1000
            file_len = len(file_text)

            #### we're looking for the end, but we take steps with length ``stepsize`` to go faster
            while match[-1 * (stepsize + len(end)):-1].count(end) == 0 and match.count("\n") < max_len:
                match = match + file_text[i:i+stepsize]
                i += stepsize

                if i > file_len:
                    break

            match = match.split(end)[0]

            # special geometry handling :/
            if idx == 3:
                if len(block_matches[3]) == len(word_matches[0]):
                    block_matches[3].append(match)
                else:
                    block_matches[3][-1] = match

            else:
                block_matches[idx].append(match)

    del file_text # here, have your RAM back!

    if len(block_matches[1]) == 0:
        raise ValueError(f"Can't find a title block - something is wrong with {filename}!")

    #### and from here, we're off to the races!
    n, g = parse_geometry(block_matches[3])
    title, link0, route_card, footer, job_types = parse_header_footer(block_matches[0], block_matches[1], block_matches[2], block_matches[4])
    energies, scf_iterations = parse_energies(word_matches[0])
    success, elapsed_time = parse_success_elapsed_time(word_matches[2], word_matches[3])
    charge, multip = parse_charge_multiplicity(word_matches[4])
    bonds = parse_bonds(block_matches[5])

    f = cctk.GaussianFile(job_types=job_types, route_card=route_card, link0=link0, footer=footer, success=success, elapsed_time=elapsed_time, title=title)

    molecules = [None] * len(g)
    properties = [{} for _ in range(len(g))]
    for idx, geom in enumerate(g):
        molecules[idx] = cctk.Molecule(n[0], geom, charge=charge, multiplicity=multip, bonds=bonds, checks=False)
        if idx < len(energies):
            properties[idx]["energy"] = energies[idx]
        if idx < len(scf_iterations):
            properties[idx]["scf_iterations"] = scf_iterations[idx]
        properties[idx]["link1_idx"] = link1idx
        properties[idx]["filename"] = filename
        properties[idx]["iteration"] = idx

    if cctk.GaussianJobType.OPT in job_types:
        rms_forces = extract_parameter(word_matches[5], 2)
        rms_disp = extract_parameter(word_matches[6], 2)

        if extended_opt_info:
            max_forces = extract_parameter(word_matches[7], 2)
            max_disp = extract_parameter(word_matches[8], 2)
            rms_grad = extract_parameter(word_matches[9], 5)
            max_grad = extract_parameter(word_matches[9], 3)
            rms_int = extract_parameter(word_matches[10], 5)
            max_int = extract_parameter(word_matches[10], 3)
            delta_e = extract_parameter(word_matches[11], 3, cast_to_float=False)

        for idx, force in enumerate(rms_forces):
            properties[idx]["rms_force"] = force
            properties[idx]["rms_displacement"] = rms_disp[idx]

            if extended_opt_info:
                if idx < len(max_forces):
                    properties[idx]["max_force"] = max_forces[idx]

                if idx < len(max_disp):
                    properties[idx]["max_displacement"] = max_disp[idx]

                if idx < len(max_grad):
                    properties[idx]["max_gradient"] = max_grad[idx]

                if idx < len(rms_grad):
                    properties[idx]["rms_gradient"] = rms_grad[idx]

                if idx < len(max_int):
                    properties[idx]["max_internal_force"] = max_int[idx]

                if idx < len(rms_int):
                    properties[idx]["rms_internal_force"] = rms_int[idx]

                if idx < len(delta_e):
                    change_in_energy = re.sub(r"Energy=", "", delta_e[idx])
                    properties[idx]["predicted_change_in_energy"] = float(change_in_energy.replace('D', 'E'))

    if cctk.GaussianJobType.FREQ in job_types:
        enthalpies = extract_parameter(word_matches[12], 6)
        if len(enthalpies) == 1:
            properties[-1]["enthalpy"] = enthalpies[0]
        elif len(enthalpies) > 1:
            raise ValueError(f"unexpected # of enthalpies found!\nenthalpies = {enthalpies}")

        gibbs_vals = extract_parameter(word_matches[13], 7)
        if len(gibbs_vals) == 1:
            properties[-1]["gibbs_free_energy"] = gibbs_vals[0]
        elif len(gibbs_vals) > 1:
            raise ValueError(f"unexpected # gibbs free energies found!\ngibbs free energies = {gibbs_vals}")

        frequencies = []
        try:
            frequencies += extract_parameter(word_matches[14], 2)
            frequencies += extract_parameter(word_matches[14], 3)
            frequencies += extract_parameter(word_matches[14], 4)
            properties[-1]["frequencies"] = sorted(frequencies)
        except Exception as e:
            raise ValueError("error finding frequencies")

        temperature = extract_parameter(word_matches[15], 1)
        if len(temperature) == 1:
            properties[-1]["temperature"] = temperature[0]
            corrected_free_energy = get_corrected_free_energy(gibbs_vals[0], frequencies, frequency_cutoff=100.0, temperature=temperature[0])
            properties[-1]["quasiharmonic_gibbs_free_energy"] = float(corrected_free_energy)

    if cctk.GaussianJobType.NMR in job_types:
        nmr_shifts = read_nmr_shifts(word_matches[16], molecules[0].num_atoms())
        if nmr_shifts is not None:
            properties[-1]["isotropic_shielding"] = nmr_shifts.view(cctk.OneIndexedArray)

        if re.search("nmr=mixed", f.route_card, flags=re.IGNORECASE) or re.search("nmr=spinspin", f.route_card,flags=re.IGNORECASE):
            couplings = read_j_couplings(block_matches[6], molecules[0].num_atoms())
            if couplings is not None:
                properties[-1]["j_couplings"] = couplings

    if cctk.GaussianJobType.FORCE in job_types:
        assert len(molecules) == 1, "force jobs should not be combined with optimizations!"
        forces = parse_forces(block_matches[7])
        properties[0]["forces"] = forces

    if cctk.GaussianJobType.POP in job_types:
        if re.search("hirshfeld", f.route_card) or re.search("cm5", f.route_card):
            charges, spins = parse_hirshfeld(block_matches[8])
            properties[-1]["hirshfeld_charges"] = charges
            properties[-1]["hirshfeld_spins"] = spins

    try:
        charges, dipole = parse_charges_dipole(block_matches[9], block_matches[10])
        properties[-1]["mulliken_charges"] = charges
        properties[-1]["dipole_moment"] = dipole
    except Exception as e:
        pass

    for mol, prop in zip(molecules, properties):
        f.ensemble.add_molecule(mol, properties=prop)

    f.check_has_properties()
    return f

def parse_geometry(blocks):
    nums = []
    geoms = []
    for block in blocks:
        current_nums = []
        current_geoms = []
        for line in block.split("\n")[4:-2]:
            if re.search("Distance", line):
                break
            pieces = list(filter(None, line.split(" ")))
            if len(pieces) != 6:
                continue
            current_nums.append(int(pieces[1]))
            current_geoms.append([float(pieces[3]), float(pieces[4]), float(pieces[5])])
        nums.append(current_nums)
        geoms.append(current_geoms)
    return nums, geoms

def parse_header_footer(route_block, title_block, footer_block, link0_block):
    link0 = dict()
    route_card = ""
    footer = None
    title = ""
    job_types = []

    title = title_block[0].split("\n")[2].strip()
    for line in route_block[0].split("\n"):
        route_card += line.lstrip()

    if len(footer_block) > 0:
        footer = "\n".join(list(footer_block[0].split("\n"))[1:])  # get rid of the first line
        footer = "\n".join([" ".join(list(filter(None, line.split(" ")))) for line in footer.split("\n")])

    for line in link0_block[0].split("\n"):
        if re.match(" \%", line):
            pieces = line[2:].split("=")
            link0[pieces[0]] = pieces[1]

    for name, member in cctk.GaussianJobType.__members__.items():
        if re.search(f" {member.value}", str(route_card), re.IGNORECASE):
            job_types.append(member)
    if cctk.GaussianJobType.SP not in job_types:
        job_types.append(cctk.GaussianJobType.SP)

    return title, link0, route_card, footer, job_types

def parse_energies(scf_done_block):
    energies = []
    iters = []

    for line in scf_done_block:
        pieces = list(filter(None, line.split(" ")))
        energies.append(float(pieces[4]))
        iters.append(int(pieces[7]))

    return energies, iters

def parse_success_elapsed_time(success_lines, time_lines):
    successes = len(success_lines)
    elapsed_time = 0
    for line in time_lines:
        fields = list(filter(None, line.split(" ")))
        elapsed_time += int(fields[2]) * 86400 + int(fields[4]) * 3600 + int(fields[6]) * 60 + float(fields[8])
    return successes, elapsed_time

def parse_charge_multiplicity(charge_line):
    fields = list(filter(None, charge_line[0].replace("=", " ").split(" ")))
    return int(fields[1]), int(fields[3])

def parse_bonds(bonding_block):
    if len(bonding_block) == 0:
        return None

    bond_array = []
    for line in bonding_block[0].split("\n"):
        if re.search(r"! R", line):
            pieces = list(filter(None, line.split(" ")))
            atoms = pieces[2].replace("R", "").replace("(", "").replace(")", "").split(",")
            try:
                bond_array.append([int(atoms[0]), int(atoms[1])])
            except Exception as e:
                raise ValueError(f"error parsing line - can't extract atoms!\n{line}\e{e}")
    return bond_array

def split_link1_to_text(filename):
    link1_blocks = []
    with open(filename, "r") as lines:
        current_text = ""
        for idx, line in enumerate(lines):
            current_text = current_text + line
            if re.search("Entering Link 1", line):
                link1_blocks.append(current_text)
                current_text = ""
        link1_blocks.append(current_text)
    return link1_blocks[1:] #### the first block is just a few lines

def extract_parameter(lines, position, cast_to_float=True):
    vals = []
    for line in lines:
        pieces = list(filter(None, line.split(" ")))
        if cast_to_float:
            vals.append(float(pieces[position]))
        else:
            vals.append(pieces[position])
    return vals

def parse_forces(force_block):
    forces = []
    for line in force_block[0].split("\n")[2:]:
        fields = re.split(" +", line)
        fields = list(filter(None, fields))

        if len(fields) == 5:
            forces.append([float(fields[2]), float(fields[3]), float(fields[4])])

    return cctk.OneIndexedArray(forces)

def parse_charges_dipole(mulliken_block, dipole_block):
    charges = []
    dipole = 0

    for line in mulliken_block[0].split("\n")[2:]:
        fields = re.split(" +", line)
        fields = list(filter(None, fields))

        if len(fields) == 3:
            charges.append(float(fields[2]))

    for line in dipole_block[0].split("\n")[1:]:
        fields = re.split(" +", line)
        fields = list(filter(None, fields))

        if len(fields) == 8:
            dipole = float(fields[7])
            break

    return cctk.OneIndexedArray(charges), dipole

def parse_hirshfeld(hirshfeld_block):
    charges = []
    spins = []
    for line in hirshfeld_block[0].split("\n")[2:]:
        fields = re.split(" +", line)
        fields = list(filter(None, fields))

        if len(fields) == 8:
            charges.append(float(fields[2]))
            spins.append(float(fields[3]))

    return cctk.OneIndexedArray(charges), cctk.OneIndexedArray(spins)

