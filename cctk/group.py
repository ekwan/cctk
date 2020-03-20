import sys, copy, re, math
import numpy as np
import networkx as nx

from abc import abstractmethod

from cctk import Molecule, OneIndexedArray
from cctk.helper_functions import get_covalent_radius, compute_angle_between, compute_rotation_matrix


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
        group = Group(attach_to, atomic_numbers=molecule.atomic_numbers, geometry=molecule.geometry, bonds=molecule.bonds.edges())
        return group

    def add_attachment_point(self, attach_to):
        """
        Adds ``attach_to``and ``adjacent`` attributes to the instance.

        Automatically centers atom ``adjacent`` on the origin, to simplify downstream mathematics.
        """
        if len(super().get_adjacent_atoms(attach_to)) != 1:
            raise ValueError(f"atom {attach_to} is making too many bonds!")

        self.attach_to = attach_to

        adjacent = super().get_adjacent_atoms(attach_to)
        assert len(adjacent) == 1, "can't substitute an atom with more than one adjacent atom!"
        self.adjacent = adjacent[0]

        adj_v = super().get_vector(self.adjacent)
        super().translate_molecule(-adj_v)

    @abstractmethod
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
        #### this code can be a bit complex: for an example, let's imagine converting benzene to toluene by adding methane (Group) to benzene (Molecule)
        ####     add_to would be the benzene H (atom on Molecule you replace with the new group)
        ####     adjacent_atom would be the benzene C
        ####     group.attach_to would be the methane H
        ####     group.adjacent would be the methane C

        #### prevent in-place modification of molecule - could lead to pernicious errors!

        try:
            add_to = int(add_to)
        except:
            raise TypeError("add_to not castable to int")

        molecule = copy.deepcopy(molecule)
        molecule._check_atom_number(add_to)

        adjacent_atom = molecule.get_adjacent_atoms(add_to)
        assert (
            len(adjacent_atom) > 0
        ), "can't substitute an atom without an adjacent atom! (are there bonds defined for this molecule? consider calling molecule.assign_connectivity()!)"
        assert len(adjacent_atom) == 1, "can't substitute an atom with more than one adjacent atom!"
        adjacent_atom = adjacent_atom[0]

        attach_to = group.attach_to
        other_indices = np.ones_like(group.atomic_numbers).astype(bool)
        other_indices[attach_to] = False 
        other_indices[group.adjacent] = False

        #### we need to change the bond length somewhat to prevent strange behavior
        old_radius = get_covalent_radius(molecule.atomic_numbers[add_to])
        new_radius = get_covalent_radius(group.atomic_numbers[group.adjacent])
        delta_rad = new_radius - old_radius

        #### make the swap! (this only adds the atoms, still have to get the geometry right)
        molecule.atomic_numbers[add_to] = group.atomic_numbers[group.adjacent]
        new_indices = [i + molecule.num_atoms() for i in range(1, np.sum(other_indices) + 1)]
        molecule.atomic_numbers = np.hstack([molecule.atomic_numbers, group.atomic_numbers[other_indices]])
        molecule.atomic_numbers = molecule.atomic_numbers.view(OneIndexedArray)

        #### have to keep track of what all the new indices are, to carry over connectivity
        new_indices.insert(group.adjacent - 1, add_to)
        new_indices.insert(attach_to - 1, adjacent_atom)

        #### adjust the bond length by moving add_to
        molecule.set_distance(adjacent_atom, add_to, molecule.get_distance(adjacent_atom, add_to) + delta_rad)

        #### rotate group to match the new positioning
        v_g = group.get_vector(group.attach_to, group.adjacent)
        v_m = molecule.get_vector(add_to, adjacent_atom)
        theta = compute_angle_between(v_g, v_m)

        #### rotate each atom and add it...
        center_pos = molecule.get_vector(add_to)
        rot = compute_rotation_matrix(np.cross(v_g, v_m), -(180 - theta))
        for vector in group.geometry[other_indices]:
            new_v = np.dot(rot, vector) + center_pos
            molecule.geometry = np.vstack((molecule.geometry, new_v))
            molecule.geometry = molecule.geometry.view(OneIndexedArray)

        #### now we have to merge the new bonds
        for (atom1, atom2) in group.bonds.edges():
            molecule.add_bond(new_indices[atom1-1], new_indices[atom2-1])

        assert len(molecule.atomic_numbers) == len(
            molecule.geometry
        ), f"molecule has {len(molecule.atomic_numbers)} atoms but {len(molecule.geometry)} geometry elements!"


        #### now we want to find the "lowest" energy conformation, defined as the rotamer which minimizes the RMS distance between all atoms
        if group.num_atoms() > 3:
            adjacent_on_old_molecule = molecule.get_adjacent_atoms(adjacent_atom)[0]
            adjacent_on_new_molecule = molecule.get_adjacent_atoms(add_to)[-1]
            molecule.optimize_dihedral(adjacent_on_old_molecule, adjacent_atom, add_to, adjacent_on_new_molecule)

        try:
            molecule.check_for_conflicts()
        except ValueError as error_msg:
            raise ValueError(f"molecule contains conflicts: {str(error_msg)}!")

        return molecule
