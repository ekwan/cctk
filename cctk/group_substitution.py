import numpy as np
import math
import re
import copy

from cctk import Molecule, Group
from cctk.helper_functions import get_covalent_radius, compute_angle_between, compute_rotation_matrix

#### Helper file to deal with group substitution, since the code gets a bit hairy

def add_group_to_molecule(molecule, group, add_to):
    """
    Adds a `Group` object to a `Molecule` at the specified atom, and returns a new `Molecule` object (generated using `copy.deepcopy()`).
    Automatically attempts to detect clashes by rotating group until no clashes are found

    The atom in `group` that replaces `add_to` in `molecule` will inherit the number of `add_to` - however, the other atoms in `group` will be appended to the atom list.

    Args:
        molecule (Molecule): the molecule to change
        group (Group): the group to affix
        add_to (int): the 1-indexed atom number on `molecule` to add `group` to
    """

    #### prevent in-place modification of molecule - could lead to pernicious errors!
    molecule = copy.deepcopy(molecule)
    molecule._check_atom_number(add_to)

    adjacent_atom = molecule.get_adjacent_atoms(add_to)
    assert (len(adjacent_atom) == 1), "can't substitute an atom with more than one adjacent atom!"
    adjacent_atom = adjacent_atom[0]

    attach_to = group.attach_to
    other_indices = np.ones_like(group.atoms)
    other_indices[attach_to-1] = 0

    #### we need to change the bond length somewhat to prevent strange behavior
    old_radius = get_covalent_radius(molecule.atomic_numbers[add_to-1])
    new_radius = get_covalent_radius(group.atomic_numbers[attach_to-1])
    delta_rad = new_radius - old_radius

    #### make the swap! (this only adds the atoms, still have to get the geometry right)
    molecule.atoms[add_to-1] = group.atomic_numbers[attach_to-1]
    new_indices = [i + molecule.num_atoms() for i in range(len(other_indices))]
    molecule.atoms.append(group.atoms[other_indices])

    #### adjust the bond length by moving add_to
    molecule.set_distance(adjacent_atom, add_to) = molecule.get_distance(adjacent_atom, add_to) + delta_rad

    #### rotate group to match the new positioning
    v_g = group.get_vector(group.attach_to, group.adjacent)
    v_m = molecule.get_vector(molecule.add_to, adjacent_atom)
    theta = compute_angle_between(v_g, v_m)

    new_center = molecule.get_vector(molecule.add_to)
    rot = group.compute_rotation_matrix(np.cross(v_g,v_m), -theta)
    for vector in group.geometry[other_indices]:
        new_v = np.dot(rot, vector) + new_center
        molecule.geometry.append(new_v)

    assert (len(molecule.atomic_numbers) == len(molecule.geometry)), f"molecule has {len(molecule.atomic_numbers)} atoms but {len(molecule.geometry)} geometry elements!"

    #### checks for conflicts, etc

    return molecule
