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

def read_file_fast(file_text, filename, link1idx, max_len=10000, extended_opt_info=False):

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
            "Leave Link  202",
            1000,
        ],
        ["Wallingford", "#p", 1],
        ["Initial Parameters", "! A", 1], #5
        ["Total nuclear spin-spin coupling J", "Leave Link", 1],
        ["Forces (Hartrees/Bohr)", "Cartesian Forces", 1],
        ["Hirshfeld charges, spin densities, dipoles, and CM5 charges", " Hirshfeld charges", 1],
        ["Mulliken charges", "Sum of Mulliken charges", 1],
        ["Electronic spatial extent", "Quadrupole moment", 1], #10
        ["normal coordinates", "Thermochemistry", 1],
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

        vibrational_modes = parse_modes(block_matches[11], num_atoms=molecules[-1].num_atoms(), hpmodes=re.search("hpmodes", route_card))
        molecules[-1].vibrational_modes = vibrational_modes

        frequencies = []
        try:
            frequencies += extract_parameter(word_matches[14], 2)

            # very small molecules might only have 1 or 2 freqs
            try:
                frequencies += extract_parameter(word_matches[14], 3)
            except:
                pass
            try:
                frequencies += extract_parameter(word_matches[14], 4)
            except:
                pass

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

def parse_modes(freq_block, num_atoms, hpmodes=False):
    freqs = list()
    masses = list()
    force_ks = list()
    displacements = list()

    if len(freq_block) == 0:
        return None

    chunks = freq_block[0].split("Freq")

    if hpmodes:
        chunks = chunks[1:]

    for chunk in chunks:
        lines = chunk.split("\n")

        if hpmodes:
            num_cols = len(re.split(" +", lines[0])) - 2
            current_displacements = [np.zeros(shape=(num_atoms, 3)) for x in range(num_cols)]

            if len(freqs):
                new_freqs = list(filter(None, re.split(" +", lines[0])))[2:]

                if float(new_freqs[-1]) <= float(freqs[-1]):
                    break # want to skip the non-hpmodes section, so no looping allowed
                else:
                    freqs += new_freqs
            else:
                freqs += list(filter(None, re.split(" +", lines[0])))[2:]

            masses += list(filter(None, re.split(" +", lines[1])))[3:]
            force_ks += list(filter(None, re.split(" +", lines[2])))[3:]

            for line in lines[6:]:
                fields = re.split(" +", line)
                fields = list(filter(None, fields))

                if len(fields) < (num_cols + 3):
                    continue

                if fields[0] == "Harmonic":
                    break

                for col_idx, val in enumerate(fields[3:]):
                    current_displacements[col_idx][int(fields[1])-1][int(fields[0])-1] = val

            for d in current_displacements:
                displacements.append(d.view(cctk.OneIndexedArray))

        else:
            current_displacements = [list() for x in re.split(" +", lines[0])[2:]]

            freqs += re.split(" +", lines[0])[2:]
            masses += re.split(" +", lines[1])[4:]
            force_ks += re.split(" +", lines[2])[4:]

            for line in lines[5:]:
                fields = re.split(" +", line)
                fields = list(filter(None, fields))

                if len(fields) < 4:
                    break

                current_displacements[0].append([float(x) for x in fields[2:5]])

                if len(current_displacements) > 1:
                    current_displacements[1].append([float(x) for x in fields[5:8]])

                if len(current_displacements) > 2:
                    current_displacements[2].append([float(x) for x in fields[8:11]])


            for d in current_displacements:
                displacements.append(cctk.OneIndexedArray(d))

    freqs = [float(x) for x in freqs]
    masses = [float(x) for x in masses]
    force_ks = [float(x) for x in force_ks]

    assert len(freqs) == len(masses)
    assert len(freqs) == len(force_ks)
    assert len(freqs) == len(displacements)

    modes = list()
    for f, m, k, d in zip(freqs, masses, force_ks, displacements):
        k *= 143.9326 # mdyne Å**-1 to kcal/mol Å**-2
        modes.append(cctk.VibrationalMode(frequency=f, reduced_mass=m, force_constant=k, displacements=d))

    return modes
