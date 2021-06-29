import copy
import numpy as np
import networkx as nx

import cctk
from cctk.helper_functions import get_covalent_radius, compute_angle_between, compute_rotation_matrix


class Group(cctk.Molecule):
    """
    Class representing a functional group.

    Note that a Group instance does not need to be missing atoms. Rather, the atom given by `attach_to` will be replaced wholesale by another molecule, and the bond distances scaled automatically.

    Attributes:
        attach_to (int): atom number to replace with larger fragment. must have only one bond! (e.g. H in F3C-H)
        adjacent (int): atom number that will be bonded to new molecule. (e.g. C in F3C-H)
        isomorphic (list of lists): list of lists of atoms that should be considered symmetry equivalent.
            For instance, the three methyl protons can be considered symmetry equivalent, so ``methane.isomorphic = [[3, 4, 5]]``.
        _map_from_truncated(dict): a dictionary mapping atom numbers of the group without ``attach_to`` to the atom numbers of the normal group
    """

    def __init__(self, attach_to, isomorphic=None, **kwargs):
        super().__init__(**kwargs)
        self.add_attachment_point(attach_to)
        self._map_from_truncated = None

        if isomorphic is not None:
            assert isinstance(isomorphic, list), "group.isomorphic must be list of lists!"
        self.isomorphic = isomorphic

    @classmethod
    def new_from_molecule(cls, molecule, attach_to, **kwargs):
        """
        Convenient method to convert ``molecule`` to ``group`` directly.
        """
        group = Group(attach_to, atomic_numbers=molecule.atomic_numbers, geometry=molecule.geometry, bonds=molecule.bonds.edges(), **kwargs)
        return group

    def add_attachment_point(self, attach_to):
        """
        Adds ``attach_to`` and ``adjacent`` attributes to the instance.

        Automatically centers atom ``adjacent`` on the origin, to simplify downstream mathematics.
        """
        n_bonds = len(super().get_adjacent_atoms(attach_to))
        if n_bonds != 1:
            raise ValueError(f"atom {attach_to} is making {n_bonds} but must make 1 bond to be a valid attachment point")

        self.attach_to = attach_to

        adjacent = super().get_adjacent_atoms(attach_to)
        assert len(adjacent) == 1, "can't substitute an atom with more than one adjacent atom!"
        self.adjacent = adjacent[0]

        adj_v = super().get_vector(self.adjacent)
        super().translate_molecule(-adj_v)

    @staticmethod
    def add_group_to_molecule(molecule, group, add_to, optimize=True, return_mapping=False):
        """
        Adds a `Group` object to a `Molecule` at the specified atom, and returns a new `Molecule` object (generated using `copy.deepcopy()`).
        Automatically attempts to prevent clashes by minimizing pairwise atomic distances.

        The atom in `group` that replaces `add_to` in `molecule` will inherit the number of `add_to` - however, the other atoms in `group` will be appended to the atom list.

        Args:
            molecule (Molecule): the molecule to change
            group (Group): the group to affix
            add_to (int): the 1-indexed atom number on `molecule` to add `group` to
            optimize (bool): whether or not to perform automated dihedral optimization
            return_mapping (bool): whether or not to return dictionaries mapping atom numbers from starting materials to products

        Returns:
            new Molecule object

            (optional) molecule_to_new dictionary mapping atom numbers from starting molecule (key) to new atom numbers (val)
            (optional) group_to_new dictionary mapping atom numbers from starting group (key) to new atom numbers (val)
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
        original_num_atoms = molecule.num_atoms()

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
        molecule.atomic_numbers = molecule.atomic_numbers.view(cctk.OneIndexedArray)

        #### have to keep track of what all the new indices are, to carry over connectivity
        new_indices.insert(group.adjacent - 1, add_to)
        new_indices.insert(attach_to - 1, adjacent_atom)

        #### track atom number mapping
        molecule_to_new = {z : z for z in range(1, molecule.num_atoms() + 1)}
        molecule_to_new[add_to] = None

        group_to_new = {}
        offset = 1
        for z in range(1, group.num_atoms() + 1):
            if other_indices[z]:
                group_to_new[z] = original_num_atoms + offset
                offset += 1
            else:
                group_to_new[z] = None
        group_to_new[group.adjacent] = add_to

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
            molecule.geometry = molecule.geometry.view(cctk.OneIndexedArray)

        #### now we have to merge the new bonds
        for (atom1, atom2) in group.bonds.edges():
            molecule.add_bond(new_indices[atom1-1], new_indices[atom2-1])
        assert molecule.get_bond_order(add_to, adjacent_atom), "we didn't add the bond we were supposed to form!"

        assert len(molecule.atomic_numbers) == len(
            molecule.geometry
        ), f"molecule has {len(molecule.atomic_numbers)} atoms but {len(molecule.geometry)} geometry elements!"

        #### now we want to find the "lowest" energy conformation, defined as the rotamer which minimizes the RMS distance between all atoms
        if group.num_atoms() > 3 and optimize:
            adjacent_on_old_molecule = molecule.get_adjacent_atoms(adjacent_atom)[0]
            adjacent_on_new_molecule = molecule.get_adjacent_atoms(add_to)[-1]
            molecule.optimize_dihedral(adjacent_on_old_molecule, adjacent_atom, add_to, adjacent_on_new_molecule)

        if molecule.check_for_conflicts():
            if return_mapping:
                return molecule, molecule_to_new, group_to_new
            else:
                return molecule
        else:
            raise ValueError(f"molecule contains conflicts!")

    @staticmethod
    def remove_group_from_molecule(molecule, atom1, atom2, return_mapping=False):
        """
        The microscopic reverse of ``add_group_to_molecule`` -- splits a ``Molecule`` along the ``atom1``–``atom2`` bond
        and returns a new ``Molecule`` object (the ``atom1`` side) and a new ``Group`` (the ``atom2`` side).

        The new objects will be capped with hydrogens; atom ordering will be preserved!

        Args:
            molecule (Molecule): the molecule to change
            atom1 (int): the 1-indexed atom number on `molecule` to make part of the new ``Molecule`` object
            atom2 (int): the 1-indexed atom number on `molecule` to make part of the new ``Group`` object
            return_mapping (bool): whether or not to return dictionaries mapping atom numbers from starting materials to products

        Returns:
            new Molecule object
            new Group object

            (optional) molecule_to_molecule dictionary mapping atom numbers from starting molecule (key) to new molecule atom numbers (val)
            (optional) molecule_to_group dictionary mapping atom numbers from starting molecule (key) to new group atom numbers (val)
        """
        try:
            atom1 = int(atom1)
            atom2 = int(atom2)
        except:
            raise TypeError("atom numbers not castable to int")

        molecule = copy.deepcopy(molecule)
        molecule._check_atom_number(atom1)
        molecule._check_atom_number(atom2)

        #### define mapping dicts
        fragment1, fragment2 = molecule._get_bond_fragments(atom1, atom2)
        molecule_to_molecule = {x: i+1 for i, x in enumerate(fragment1)}
        molecule_to_group = {x: i+1 for i, x in enumerate(fragment2)}

        #### create new molecules
        new_mol = cctk.Molecule(molecule.atomic_numbers[fragment1], molecule.geometry[fragment1])
        group = cctk.Molecule(molecule.atomic_numbers[fragment2], molecule.geometry[fragment2])

        #### add capping H to new_mol
        new_mol.add_atom("H", molecule.geometry[atom2])
        molecule_to_molecule[atom2] = new_mol.num_atoms()
        old_radius = get_covalent_radius(molecule.atomic_numbers[atom2])
        H_radius = get_covalent_radius(1)
        new_dist = new_mol.get_distance(molecule_to_molecule[atom1], molecule_to_molecule[atom2]) - old_radius + H_radius
        new_mol.set_distance(molecule_to_molecule[atom1], molecule_to_molecule[atom2], new_dist)
        new_mol.add_bond(molecule_to_molecule[atom1], molecule_to_molecule[atom2])

        #### add capping H to new group
        group.add_atom("H", molecule.geometry[atom1])
        molecule_to_group[atom1] = group.num_atoms()
        old_radius = get_covalent_radius(molecule.atomic_numbers[atom1])
        new_dist = group.get_distance(molecule_to_group[atom2], molecule_to_group[atom1]) - old_radius + H_radius
        group.set_distance(molecule_to_group[atom2], molecule_to_group[atom1], new_dist)
        group.add_bond(molecule_to_group[atom2], molecule_to_group[atom1])

        #### add bonds to nascent molecules
        molecule.remove_bond(atom1, atom2)
        for (a1, a2) in molecule.bonds.edges():
            if a1 in fragment1:
                assert a2 in fragment1, "somehow we have another bond between the two groups!"
                assert molecule_to_molecule[a1] is not None, f"we don't have a mapping for atom {a1}"
                assert molecule_to_molecule[a2] is not None, f"we don't have a mapping for atom {a2}"
                new_mol.add_bond(molecule_to_molecule[a1], molecule_to_molecule[a2])
            elif a2 in fragment2:
                assert a2 in fragment2, "somehow we have another bond between the two groups!"
                assert molecule_to_group[a1] is not None, f"we don't have a mapping for atom {a1}"
                assert molecule_to_group[a2] is not None, f"we don't have a mapping for atom {a2}"
                group.add_bond(molecule_to_group[a1], molecule_to_group[a2])

        #### create Group object from group
        group = cctk.Group.new_from_molecule(attach_to=molecule_to_group[atom1], molecule=group)

        if return_mapping:
            return new_mol, group, molecule_to_molecule, molecule_to_group
        else:
            return new_mol, group

    def map_from_truncated(self):
        """
        Returns a dictionary mapping atomic numbers without ``attach_to`` to atomic_numbers with ``attach_to``.
        """
        if self._map_from_truncated is not None:
            return self._map_from_truncated

        assert self.bonds.number_of_edges() > 0, "need a bond graph to perform this operation -- try calling self.assign_connectivity()!"
        g = copy.deepcopy(self)
        g._add_atomic_numbers_to_nodes()
        tg = copy.deepcopy(g)
        tg.remove_atom(g.attach_to)

        nm = nx.algorithms.isomorphism.categorical_node_match("atomic_number", 0)
        match = nx.algorithms.isomorphism.GraphMatcher(g.bonds, tg.bonds, node_match=nm)

        for sg in match.subgraph_isomorphisms_iter():
            if self.attach_to in sg.keys():
                continue
            sg = {v: k for k, v in sg.items()} # invert
            self._map_from_truncated = sg
            return sg
