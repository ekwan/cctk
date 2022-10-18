import math, copy, re
import numpy as np
import networkx as nx
from scipy.spatial.distance import cdist
import pkg_resources
import yaml

import cctk
from cctk.helper_functions import (
    get_symbol,
    get_number,
    compute_rotation_matrix,
    compute_distance_between,
    compute_angle_between,
    compute_dihedral_between,
    compute_unit_vector,
    get_covalent_radius,
    get_vdw_radius,
    numpy_to_bytes,
    bytes_to_numpy,
    _recurse_through_formula,
)
import cctk.topology as top

class Molecule:
    """
    Class representing a single molecular geometry.

    In contrast to typical Python behavior, ``atomic_numbers`` and ``geometry`` are indexed from one, to simplify interfacing with computational chemistry programs.
    This has been done by defining a custom wrapper for ``numpy.ndarray`` called ``cctk.OneIndexedArray``.

    All other datatypes are indexed from 0.

    Attributes:
        name (str): for identification, optional
        atomic_numbers (cctk.OneIndexedArray, dtype=np.int8): list of atomic numbers
        geometry (cctk.OneIndexedArray): list of 3-tuples of xyz coordinates - same ordering as ``atomic_numbers``
        bonds (nx.Graph or list of tuples): connectivity graph or list of 2-tuples, with each element representing the 1-indexed atom number of a bonded pair
        charge (int): the charge of the molecule
        multiplicity (int): the spin state of the molecule (1 corresponds to singlet, 2 to doublet, 3 to triplet, etc. -- so a multiplicity of 1 is equivalent to S=0)
        vibrational_modes (list of cctk.VibrationalMode): vibrational modes
    """

    def __init__(self, atomic_numbers, geometry, name=None, bonds=None, charge=0, multiplicity=1, checks=True):
        """
        Create new Molecule object, and assign connectivity if needed.

        ``bonds`` must be a list of edges (i.e. an n x 2 ``numpy`` array).

        If ``checks`` is True, the atomic numbers in bonds will all be checked for consistency.
        This option can be disabled by setting ``checks`` to False, but this is not recommended for external data.
        """
        if len(atomic_numbers) != len(geometry):
            raise ValueError(f"length of geometry ({len(geometry)}) and atomic_numbers ({len(atomic_numbers)}) does not match!\n{atomic_numbers}\n{geometry}")

        try:
            atomic_numbers = np.asarray(atomic_numbers, dtype=np.int8).view(cctk.OneIndexedArray)
        except Exception as e:
            raise ValueError("invalid atom list")

        try:
            geometry = np.array(geometry, dtype=np.float32).view(cctk.OneIndexedArray)
        except Exception as e:
            raise TypeError("geometry cannot be cast to ``np.ndarray`` of floats!")

        if name is not None:
            if not isinstance(name, str):
                raise TypeError("name must be a string!")

        if not isinstance(charge, int):
            try:
                charge = int(charge)
            except Exception as e:
                raise TypeError("charge must be integer or castable to integer!")

        if not isinstance(multiplicity, int):
            try:
                multiplicity = int(multiplicity)
            except Exception as e:
                raise TypeError("multiplicity must be positive integer or castable to positive integer")
        assert multiplicity > 0, "multiplicity must be positive"

        self.atomic_numbers = atomic_numbers
        self.geometry = geometry

        self.name = name
        self.multiplicity = multiplicity
        self.charge = charge

        self.vibrational_modes = list()

        if isinstance(bonds, nx.Graph):
            self.bonds = bonds
        elif isinstance(bonds, (list,np.ndarray,nx.classes.reportviews.EdgeView)):
            if checks:
                known_atomic_numbers = set()
                for bond in bonds:
                    assert len(bond)==2, f"unexpected number of atoms in bond, expected 2, got {len(bond)}"
                    if bond[0] not in known_atomic_numbers:
                        self._check_atom_number(bond[0])
                        known_atomic_numbers.add(bond[0])
                    if bond[1] not in known_atomic_numbers:
                        self._check_atom_number(bond[1])
                        known_atomic_numbers.add(bond[1])

            self.bonds = nx.Graph()
            self.bonds.add_nodes_from(range(1, len(atomic_numbers) + 1))
            self.bonds.add_edges_from(bonds, weight=1)
        elif bonds is None:
            self.bonds = nx.Graph()
            self.bonds.add_nodes_from(range(1, len(atomic_numbers)+1))
        else:
            raise ValueError(f"unexpected type for bonds: {type(bonds)}")

    def __str__(self):
        if self.name is not None:
            return f"Molecule (name={self.name}, {len(self.atomic_numbers)} atoms)"
        else:
            return f"Molecule ({len(self.atomic_numbers)} atoms)"

    def __repr__(self):
        return str(self) # placeholder

#    def __eq__(self, other):
    @classmethod
    def equal(cls, mol1, mol2):
        """
        Atomic numbers, geometry, charge, and multiplicity all must match. Name is irrelevant.
        """
        if not isinstance(mol1, cctk.Molecule):
            return False

        if not isinstance(mol2, cctk.Molecule):
            return False

        comparisons = [
            np.array_equal(mol1.atomic_numbers, mol2.atomic_numbers),
            np.array_equal(mol1.geometry, mol2.geometry),
            mol1.charge == mol2.charge,
            mol1.multiplicity == mol2.multiplicity
        ]

        return all(comparisons)

    def assign_connectivity(self, cutoff=0.2, periodic_boundary_conditions=None):
        """
        Automatically recalculates bonds based on covalent radii. If two atoms are closer than the sum of their covalent radii + ``cutoff`` Angstroms,
        then they are considered bonded.

        Args:
            cutoff (float): the threshold (in Angstroms) for how close two covalent radii must be to be considered bonded

        Returns:
            self
        """

        #### delete all edges
        self.bonds = nx.create_empty_copy(self.bonds)

        assert isinstance(cutoff, (float, int)), "need cutoff to be numeric!"
        g = self.geometry.view(np.ndarray)

        dist_matrix = None

        #### cdist is SO FAST
        if periodic_boundary_conditions is None:
            dist_matrix = cdist(g, g, "euclidean")
        else:
            # even 16 cdist calls is faster than any other implementation, i tested it
            pbc = periodic_boundary_conditions
            assert isinstance(pbc, np.ndarray) and len(pbc) == 3, "Need 3-element ``np.ndarray`` for PBCs"

            nearby_cells = [
                [0, 0, 0],
                [pbc[0], 0, 0],
                [0, pbc[1], 0],
                [0, 0, pbc[2]],
                [pbc[0], pbc[1], 0],
                [pbc[0], 0, pbc[2]],
                [0, pbc[1], pbc[2]],
                [pbc[0], pbc[1], pbc[2]],
            ]

            dist_matrices = [cdist(g, g + np.array(nc), "euclidean") for nc in nearby_cells]
            dist_matrices += [cdist(g, g - np.array(nc), "euclidean") for nc in nearby_cells]
            distances_3d = np.stack(dist_matrices)
            dist_matrix = distances_3d.min(axis=0)

        covalent_radii = {z: get_covalent_radius(z) for z in set(self.atomic_numbers)}
        radii_by_num = [covalent_radii[z] for z in self.atomic_numbers]

        for i in range(1, self.num_atoms() + 1):
            r_i = radii_by_num[i-1]
            for j in range(i + 1, self.num_atoms() + 1):
                distance = dist_matrix[i-1][j-1]
                r_j = radii_by_num[j-1]

                # 0.5 A distance is used by RasMol and Chime (documentation available online) and works well, empirically
                if distance < (r_i + r_j + cutoff):
                    self.add_bond(i, j)

        return self

    def check_for_conflicts(self, min_buffer=1, group1=None, group2=None):
        """
        Automatically checks for conflicts based on covalent radii. If two atoms are closer than the sum of their covalent radii + buffer, then they are considered clashing.
        If `group1` and `group2` are selected, then conflicts will only be evaluated between these two groups of atoms.

        Args:
            min_buffer (float): the threshold (in Angstroms) for how close two covalent radii must be to be considered clashing. 1.0 A is default, empirically.
            group1 (list): atoms to evaluate against `group2` (if `None`, defaults to all atoms)
            group2 (list): atoms to evaluate against `group1` (if `None`, defaults to all atoms)

        Returns:
            True if there are no conflicts
            ValueError if there is a conflict
        """

        if group1 is None:
            group1 = list(range(1, self.num_atoms() + 1))

        if group2 is None:
            group2 = list(range(1, self.num_atoms() + 1))

        for atom in group1 + group2:
            self._check_atom_number(atom)

        for i in group1:
            for j in group2:
                if i == j:
                    continue
                distance = self.get_distance(i, j, check=False)
                r_i = get_covalent_radius(self.get_atomic_number(i))
                r_j = get_covalent_radius(self.get_atomic_number(j))

                # 0.5 A distance is used by RasMol and Chime (documentation available online) and works well, empirically
                if distance < (r_i + r_j - min_buffer):
#                    raise ValueError(f"atoms {i} and {j} are too close - distance {distance} A!")
                    return False

        return True

    def add_bond(self, atom1, atom2, bond_order=1, check=True):
        """
        Adds a new bond to the bond graph, or updates the existing bond order. Will not throw an error if the bond already exists.

        Args:
            atom1 (int): the number of the first atom
            atom2 (int): the number of the second atom
            bond_order (int): bond order of bond between atom1 and atom2
        """
        if check:
            self._check_atom_number(atom1)
            self._check_atom_number(atom2)
            assert isinstance(bond_order, int), f"bond order {bond_order} must be an integer"
            assert bond_order >= 0, f"bond order {bond_order} must be positive"

        if self.bonds.has_edge(atom1, atom2):
            if bond_order == 0:
                self.bonds.remove_edge(atom1, atom2)
            else:
                if self.bonds[atom1][atom2]["weight"] != bond_order:
                    self.bonds[atom1][atom2]["weight"] = bond_order
        elif bond_order > 0:
            self.bonds.add_edge(atom1, atom2, weight=bond_order)

    def remove_bond(self, atom1, atom2):
        """
        Alias for ``self.add_bond(atom1, atom2, bond_order=0)`` -- more intuitive nomenclature.
        """
        self.add_bond(atom1, atom2, bond_order=0)

    def _check_atom_number(self, number):
        """
        Helper method which performs quick checks on the validity of a given atom number.
        """
        assert isinstance(number, int), "atomic number must be integer"
        assert 0 < number <= self.num_atoms(), "atom number {number} too large! (or too small - needs to be >0)"

    def formula(self, return_dict=False):
        """
        Returns the atomic formula.

        If ``return_dict`` is ``True``, then returns a ``dictionary`` with keys elemental symbols and values the number of occurrences.

        For instance, ``water.formula()`` would return ``{'O': 1, 'H': 2}``.

        If ``return_dict`` is ``False``, then returns a stringified version of the formula according to standard rules.

        For instance, ``water.formula()`` would return ``H2O``.

        Args:
            return_dict (Bool): if the method should return a string or a dictionary

        Returns:
            a dictionary or string representing the molecule's formula
        """

        formula_dict = {}
        for atom in self.atomic_numbers:
            symbol = get_symbol(atom)
            if symbol in formula_dict:
                formula_dict[symbol] += 1
            else:
                formula_dict[symbol] = 1
        if return_dict == True:
            return formula_dict
        else:
            formula = ""
            elements = list(formula_dict.keys())

            #### H and C always come first
            if "C" in elements:
                elements.remove("C")
                formula += f"C{formula_dict['C']}"

            if "H" in elements:
                elements.remove("H")
                formula += f"H{formula_dict['H']}"

            for element in sorted(elements):
                formula += f"{element}{formula_dict[element]}"

            return formula

    def _get_bond_fragments(self, atom1, atom2):
        """
        Returns the pieces of a molecule that one would obtain by ereaking the bond between two atoms. Will throw ``ValueError`` if the atoms are in a ring.
        Useful for distance/angle/dihedral scans -- one fragment can be moved and the other held constant.

        Args:
            atom1 (int): the number of the first atom
            atom2 (int): the number of the second atom

        Returns:
            fragment1: the list of atoms in fragment 1 (containing atom1)
            fragment2: the list of atoms in fragment 2 (containing atom2)

        """

        self._check_atom_number(atom1)
        self._check_atom_number(atom2)

        assert self.bonds.number_of_edges() > 0, "need a bond graph to perform this operation -- try calling self.assign_connectivity()!"

        bond_order = self.get_bond_order(atom1, atom2)
        if self.bonds.has_edge(atom1, atom2):
            self.bonds.remove_edge(atom1, atom2)

            fragments = nx.connected_components(self.bonds)
            fragment1 = []
            fragment2 = []

            for fragment in fragments:
                if atom1 in fragment:
                    if atom2 in fragment:
                        self.add_bond(atom1, atom2, bond_order) # not adding back this bond causes some pretty pernicious errors
                        raise ValueError(f"Atom {atom1} and atom {atom2} are in a ring or otherwise connected!")
                    else:
                        fragment1 = fragment
                if atom2 in fragment:
                    fragment2 = fragment

            self.add_bond(atom1, atom2, bond_order)
            return list(fragment1), list(fragment2)
        else:
            raise ValueError(f"No bond between atom {atom1} and atom {atom2}!")

    def _get_fragment_containing(self, atom):
        """
        Get the fragment containing the atom with number ``atom``.

        Args:
            atom (int): the number of the atom

        Returns:
            a list of all the atoms in the fragment
        """

        self._check_atom_number(atom)

        fragments = nx.connected_components(self.bonds)

        for fragment in fragments:
            if atom in fragment:
                return list(fragment)

    def set_distance(self, atom1=None, atom2=None, distance=None, move="group", atoms=None):
        """
        Adjusts the ``atom1`` -- ``atom2`` bond length to be a fixed distance by moving atom2.

        If ``move`` is set to "group", then all atoms bonded to ``atom2`` will also be moved.

        If ``move`` is set to "atom", then only atom2 will be moved.

        Args:
            atom1 (int): the number of the first atom
            atom2 (int): the number of the second atom
            distance (float): distance in Angstroms of the final bond
            move (str): determines how fragment moving is handled
            atoms (list): 2-element list of atom numbers

        Returns:
            the Molecule object
        """

        if (atom1 is None) and (atom2 is None):
            assert isinstance(atoms, (list, np.ndarray)), "atom numbers need to come from fields or list!"
            assert len(atoms) == 2, "need 2 atom numbers to set distance"
            atom1 = atoms[0]
            atom2 = atoms[1]

        assert isinstance(distance, (float, int, np.number)), "need distance to set distance"

        self._check_atom_number(atom1)
        self._check_atom_number(atom2)

        if (not isinstance(distance, float)) or (distance < 0):
            raise ValueError(f"invalid value {distance} for distance!")

        atoms_to_move = []
        if move == "group":
            if self.get_bond_order(atom1, atom2):
                _, atoms_to_move = self._get_bond_fragments(atom1, atom2)
            else:
                atoms_to_move = self._get_fragment_containing(atom2)
        elif move == "atom":
            atoms_to_move = [atom2]
        else:
            raise ValueError(f"Invalid option {move} for parameter 'move'!")

        if (atom1 in atoms_to_move and atom2 in atoms_to_move) and move == "group":
            raise ValueError('both our atoms are connected which will preclude any movement with ``move`` set to "group"')

        current_distance = self.get_distance(atom1, atom2)

        v1 = self.get_vector(atom1)
        v2 = self.get_vector(atom2)
        vb = v2 - v1

        if np.linalg.norm(vb) - current_distance > 0.00001:
            raise ValueError(f"Error calculating bond distance!")

        #### move all the atoms
        delta = distance - current_distance
        unitv = compute_unit_vector(vb)
        for atom in atoms_to_move:
            self.geometry[atom] = self.geometry[atom] + (delta * unitv)

        #### check everything worked okay...
        v1f = self.get_vector(atom1)
        v2f = self.get_vector(atom2)
        vbf = v2f - v1f

        if np.linalg.norm(vbf) - distance > 0.001:
            new_dist = np.linalg.norm(vbf)
            raise ValueError(f"Error moving bonds -- new distance is {new_dist:.3f}. Operation failed!")

        return self

    def set_angle(self, atom1=None, atom2=None, atom3=None, angle=None, move="group", atoms=None):
        """
        Adjusts the ``atom1`` -- ``atom2`` -- ``atom3`` bond angle to be a fixed value by moving ``atom3``.

        If `move` is set to "group", then all atoms bonded to ``atom3`` will also be moved.

        If `move` is set to "atom", then only ``atom3`` will be moved.

        Args:
            atom1 (int): the number of the first atom
            atom2 (int): the number of the second atom
            atom3 (int): the number of the third atom
            angle (float): final value in degrees of the ``atom1`` -- ``atom2`` -- ``atom3`` angle
            move (str): determines how fragment moving is handled
            atoms (list): 3-element list of atom numbers

        Returns:
            the Molecule object
        """

        if (atom1 is None) and (atom2 is None) and (atom3 is None) :
            assert isinstance(atoms, (list, np.ndarray)), "atom numbers need to come from fields or list!"
            assert len(atoms) == 3, "need 3 atom numbers to set angle"
            atom1 = atoms[0]
            atom2 = atoms[1]
            atom3 = atoms[2]

        assert isinstance(angle, (float, int, np.number)), "need angle to set angle"

        self._check_atom_number(atom1)
        self._check_atom_number(atom2)
        self._check_atom_number(atom3)

        if self.get_distance(atom1, atom2) < 0.01:
            raise ValueError(f"atom {atom1} and atom {atom2} are too close!")

        if self.get_distance(atom2, atom3) < 0.01:
            raise ValueError(f"atom {atom2} and atom {atom3} are too close!")

        if self.get_distance(atom1, atom3) < 0.01:
            raise ValueError(f"atom {atom1} and atom {atom3} are too close!")

        try:
            angle = float(angle)
        except Exception as e:
            raise TypeError(f"angle {angle} cannot be converted to float!")

        if (not isinstance(angle, float)) or ((angle < 0) or (angle > 360)):
            raise ValueError(f"invalid value {angle} for angle!")

        atoms_to_move = []
        if move == "group":
            if self.get_bond_order(atom2, atom3):
                _, atoms_to_move = self._get_bond_fragments(atom2, atom3)
            elif self.are_connected(atom2, atom3):
                raise ValueError(
                    f"atom {atom2} and atom {atom3} are connected but not bonded -- cannot adjust angle! try manually removing one or more bonds."
                )
            else:
                atoms_to_move = self._get_fragment_containing(atom3)
        elif move == "atom":
            atoms_to_move = [atom3]
        else:
            raise ValueError(f"Invalid option {move} for parameter 'move'!")

        if atom1 in atoms_to_move:
            raise ValueError(
                f"atom {atom1} and atom {atom3} are connected in multiple ways -- cannot adjust angle! try manually removing one or more bonds."
            )

        current_angle = self.get_angle(atom1, atom2, atom3)
        delta = angle - current_angle

        if np.abs(delta) < 0.001:
            return

        #### now the real work begins...

        #### move everything to place atom2 at the origin
        v2 = self.get_vector(atom2)
        self.translate_molecule(-v2)

        v1 = self.get_vector(atom1)
        v3 = self.get_vector(atom3)

        #### perform the actual rotation
        rot_axis = np.cross(v1, v3)
        rot_matrix = compute_rotation_matrix(rot_axis, delta)
        for atom in atoms_to_move:
            self.geometry[atom] = np.dot(rot_matrix, self.get_vector(atom))

        #### and move it back!
        self.translate_molecule(v2)

        final_angle = self.get_angle(atom1, atom2, atom3)

        #### need to compare cosines to prevent insidious phase difficulties (like 0.00 and 359.99)
        if np.abs(math.cos(math.radians(final_angle)) - math.cos(math.radians(angle))) > 0.001:
            raise ValueError(f"Error rotating atoms -- expected angle {angle}, got {final_angle}  -- operation failed!")

        return self

    def set_dihedral(self, atom1=None, atom2=None, atom3=None, atom4=None, dihedral=None, move="group34", check_result=True, atoms=None):
        """
        Adjusts the ``atom1`` -- ``atom2`` -- ``atom3`` -- ``atom4`` dihedral angle to be a fixed value by moving atom 4.

        If ``move`` is set to "atom", then only ``atom4`` will be moved.

        If ``move`` is set to "group4", then all atoms bonded to ``atom4`` will also be moved.

        If ``move`` is set to "group34", then all atoms bonded to ``atom3`` and ``atom4`` will also be moved.

        Args:
            atom1 (int): the number of the first atom
            atom2 (int): the number of the second atom
            atom3 (int): the number of the third atom
            atom4 (int): the number of the fourth atom
            dihedral (float): final value in degrees of the ``atom1`` -- ``atom2`` -- ``atom3`` -- ``atom4`` angle
            move (str): determines how fragment moving is handled
            check_result (Bool): whether the final answer should be checked for correctness
            atoms (list): 4-element list of atomic numbers

        Returns:
            the Molecule object
        """

        if (atom1 is None) and (atom2 is None) and (atom3 is None) and (atom4 is None):
            assert isinstance(atoms, (list, np.ndarray)), "atom numbers need to come from fields or list!"
            assert len(atoms) == 4, "need 4 atom numbers to set dihedral"
            atom1 = atoms[0]
            atom2 = atoms[1]
            atom3 = atoms[2]
            atom4 = atoms[3]

        assert isinstance(dihedral, (float, int, np.number)), "need angle to set dihedral angle"

        # check atom numbers
        self._check_atom_number(atom1)
        self._check_atom_number(atom2)
        self._check_atom_number(atom3)
        self._check_atom_number(atom4)

        # check there is bond connectivity information
        assert len(self.bonds) > 0, "no bond connectivity information"

        # check for collinearity
        angle = self.get_angle(atom1, atom2, atom3, check=False)
        assert 0.0001 < angle < 179.9999, f"1/2/3 atoms {atom1}-{atom2}-{atom3} are collinear (angle={angle:.8f})"
        angle = self.get_angle(atom2, atom3, atom4, check=False)
        assert 0.0001 < angle < 179.9999, f"2/3/4 atoms {atom2}-{atom3}-{atom4} are collinear (angle={angle:.8f})"

        for x in [atom1, atom2, atom3, atom4]:
            for y in [atom1, atom2, atom3, atom4]:
                if x <= y:
                    continue
                else:
                    if self.get_sq_distance(x, y, check=False) < 0.001:
                        raise ValueError(f"atom {x} and atom {y} are too close!")

        try:
            dihedral = float(dihedral)
        except Exception as e:
            raise TypeError(f"dihedral angle {dihedral} cannot be converted to float!")

        if (not isinstance(dihedral, float)) or ((dihedral < 0) or (dihedral > 360)):
            raise ValueError(f"invalid value {dihedral} for dihedral angle!")

        atoms_to_move = []
        if move == "group34":
            #### add atom3's fragment to atom4
            if self.get_bond_order(atom2, atom3):
                _, atoms_to_move = self._get_bond_fragments(atom2, atom3)
            elif self.are_connected(atom2, atom3):
                raise ValueError(
                    f"atom {atom2} and atom {atom3} are connected but not bonded -- cannot adjust dihedral angle! try manually removing one or more bonds."
                )
            else:
                atoms_to_move = self._get_fragment_containing(atom3)

            #### and make sure atom4 is in there too!
            if atom4 not in atoms_to_move:
                atoms_to_move += self._get_fragment_containing(atom4)
        elif move == "group4":
            if self.get_bond_order(atom3, atom4):
                _, atoms_to_move = self._get_bond_fragments(atom3, atom4)
            elif self.are_connected(atom3, atom4):
                raise ValueError(
                    f"atom {atom3} and atom {atom4} are connected but not bonded -- cannot adjust dihedral angle! try manually removing one or more bonds."
                )
            else:
                atoms_to_move = self._get_fragment_containing(atom4)
        elif move == "atom":
            atoms_to_move = [atom4]
        else:
            raise ValueError(f"Invalid option {move} for parameter 'move'!")

        if atom1 in atoms_to_move:
            raise ValueError(
                f"atom {atom1} and atom {atom4} are connected in multiple ways -- cannot adjust dihedral angle! try manually removing one or more bonds."
            )

        if atom2 in atoms_to_move:
            raise ValueError(
                f"atom {atom2} and atom {atom4} are connected in multiple ways -- cannot adjust dihedral angle! try manually removing one or more bonds."
            )

        if atom4 not in atoms_to_move:
            raise ValueError(f"atom {atom4} is not going to be moved... this operation is doomed to fail!")

        current_dihedral = self.get_dihedral(atom1, atom2, atom3, atom4, check=False)
        delta = (dihedral - current_dihedral) % 360

        if np.abs(delta) < 0.001:
            return self

        #### now the real work begins...
        #### move everything to place atom2 at the origin
        v3 = self.get_vector(atom3, check=False)
        self.translate_molecule(-v3)

        #### perform the actual rotation
        rot_matrix = compute_rotation_matrix(-self.get_vector(atom2, check=False), delta)

        for atom in atoms_to_move:
            self.geometry[atom] = np.dot(rot_matrix, self.get_vector(atom, check=False))

        #### and move it back!
        self.translate_molecule(v3)

        if check_result:
            final_dihedral = self.get_dihedral(atom1, atom2, atom3, atom4, check=False)

            #### need to compare cosines to prevent insidious phase difficulties (like 0.00 and 359.99)
            #### this will throw ValueError for differences of about 2 degrees
            if np.abs(math.cos(math.radians(final_dihedral)) - math.cos(math.radians(dihedral))) > 0.001:
                raise ValueError(f"Error rotating atoms -- expected dihedral angle {dihedral}, got {final_dihedral}  -- operation failed!")

        return self

    def translate_molecule(self, vector):
        """
        Translates the whole molecule by the given vector.

        Args:
            vector (vector): the vector to translate by

        Returns:
            the Molecule object
        """
#        for atom in range(1, self.num_atoms() + 1):
#            self.geometry[atom] = self.geometry[atom] + vector

        self.geometry += vector

        return self

    def rotate_molecule(self, axis, degrees):
        """
        Rotates the whole molecule around the given axis.

        Args:
            axis (vector): the vector to rotate about
            theta (float): how much to rotate (in degrees)

        Returns:
            the Molecule object
        """
        rot_matrix = compute_rotation_matrix(axis, degrees)

        for atom in range(1, self.num_atoms() + 1):
            self.geometry[atom] = np.dot(rot_matrix, self.geometry[atom])

        return self

    def calculate_mass_spectrum(self, **kwargs):
        """
        Generates list of m/z values.

        Final weights rounded to one decimal point (because of low-res MS).
        """
        form_vec = np.zeros(shape=92, dtype=np.int8)
        for z in self.atomic_numbers:
            form_vec[z] += 1

        masses, weights = _recurse_through_formula(form_vec, [0], [1], **kwargs)

        new_masses, indices = np.unique(np.round(masses, decimals=1), return_inverse=True)
        new_weights = np.zeros_like(new_masses)
        for k in range(len(new_weights)):
            new_weights[k] = np.sum(weights[np.nonzero(indices == k)])
        new_weights = new_weights / np.max(new_weights)

        return new_masses, new_weights

    def add_atom_at_centroid(self, symbol, atom_numbers, weighted=False):
        """
        Adds atom with symbol ``symbol`` at the centroid of the atoms in ``atom_numbers``.

        If ``weighted`` is ``True``, then the centroid calculation will take into account the atomic numbers of the atoms in question (placing the atom closer to more massive atoms).

        Otherwise, the average is unweighted.

        Args:
            symbol (str): the atomic symbol of the atom to be added
            atom_numbers (list): which atoms to put the new atom between
            weighted (Bool): if the centroid calculation should be weighted (see above)

        Returns:
            the Molecule object
        """

        if (not isinstance(atom_numbers, list)) or (len(atom_numbers) < 2):
            raise TypeError("atom_numbers must be list with at least two elements")

        if not isinstance(symbol, str):
            raise TypeError(f"symbol {symbol} must be a string!")

        coords = [None] * len(atom_numbers)
        weights = [1] * len(atom_numbers)
        for index, atom in enumerate(atom_numbers):
            self._check_atom_number(atom)
            coords[index] = self.get_vector(atom)
            if weighted == True:
                weights[index] = self.atomic_numbers[atom]

        new_coord = list(np.average(coords, weights=weights, axis=0))
        return self.add_atom(coordinates=new_coord, symbol=symbol)

    def add_atom(self, symbol, coordinates):
        """
        Add an atom with symbol ``symbol`` at position ``coordinates``.

        Args:
            symbol (str): symbol of the atom (e.g. "Cl", "Ar", "C")
            coordinates (list): the coordinates to add

        Returns:
            the Molecule object
        """

        if (not isinstance(coordinates, (list, np.ndarray)) or (len(coordinates) != 3)):
            raise TypeError("coordinates must be list with three elements")

        if not isinstance(symbol, str):
            raise TypeError(f"symbol {symbol} must be a string!")

        number = get_number(symbol)
        self.atomic_numbers = np.append(self.atomic_numbers, [number]).astype(np.int8).view(cctk.OneIndexedArray)
        self.geometry = np.append(self.geometry, [coordinates], axis=0).view(cctk.OneIndexedArray)
        self.bonds.add_node(self.num_atoms())

        return self

    def remove_atom(self, number):
        """
        Remove the atom with number ``number``.

        Args:
            number (int): number of the atom

        Returns:
            the Molecule object
        """

        self._check_atom_number(number)

        try:
            self.bonds.remove_node(number)
            self.geometry = np.delete(self.geometry, number - 1, axis=0).view(cctk.OneIndexedArray)
            self.atomic_numbers = np.delete(self.atomic_numbers, number - 1).view(cctk.OneIndexedArray)

            #### need to renumber to fill gaps
            self.bonds = nx.convert_node_labels_to_integers(self.bonds, first_label=1, ordering="sorted")

            return self
        except Exception as e:
            raise ValueError("removing atom {number} failed!")

    def get_atomic_number(self, atom):
        """
        Get the atomic number for a given atom.

        Args:
            atom (int): number of the first atom

        Returns:
            atomic_number (int): the atomic number of that atom
        """
        self._check_atom_number(atom)
        return self.atomic_numbers[atom]

    def get_atomic_symbol(self, atom):
        """
        Get the atomic symbol for a given atom.

        Args:
            atom (int): number of the first atom

        Returns:
            atomic_symbol (str): the atomic symbol of that atom
         """
        atomic_number = self.get_atomic_number(atom)
        return get_symbol(atomic_number)

    def get_atomic_symbols(self):
        """
        Get a list of atomic symbols for this Molecule.

        Returns:
            atomic_symbols (cctk.OneIndexedArray): the atomic symbols
        """
        n_atoms = self.get_n_atoms()
        l = [ self.get_atomic_symbol(i) for i in range(1,n_atoms+1) ]
        return cctk.OneIndexedArray(l)

    def get_n_atoms(self):
        """
        Determine how many atoms are in this Molecule.

        Returns
            n_atoms (int): the number of atoms
        """
        return len(self.atomic_numbers)

    def get_vector(self, atom, atom2=None, check=True):
        """
        Get the geometry vector for a given atom. If two atoms are specified, gives the vector connecting them (from ``atom2`` to ``atom``).
        ``mol.get_vector(atom)`` is thus equivalent to ``mol.get_vector(atom, origin)``.

        Args:
            atom1 (int): number of the first atom
            atom2 (int): number of the second atom (optional)
            check (Bool): whether to validate input data (can be overridden to prevent slow double-checking)

        Returns:
            a Numpy array
        """
        if check:
            self._check_atom_number(atom)

        if atom2:
            if check:
                self._check_atom_number(atom2)
            return (self.geometry[atom] - self.geometry[atom2]).view(np.ndarray)
        else:
            return self.geometry[atom].view(np.ndarray)

    def get_distance(self, atom1=None, atom2=None, check=True, _dist=compute_distance_between, atoms=None):
        """
        Wrapper to compute distance between two atoms.

        This function is relatively slow (rate-limiting for certain applications), so performance boosts have been implemented (e.g. preloading ``_dist``).

        Args:
            atom1 (int): number of the first atom
            atom2 (int): number of the second atom
            check (Bool): whether to validate input data (can be overridden to prevent slow double-checking)
            _dist (function): function usd to compute distance
            atoms (list): list of atomic numbers

        Returns:
            the distance, in Angstroms
        """
        if (atom1 is None) and (atom2 is None):
            assert isinstance(atoms, (list, np.ndarray)), "atom numbers need to come from fields or list!"
            assert len(atoms) == 2, "need 2 atom numbers to get distance"
            atom1 = atoms[0]
            atom2 = atoms[1]

        if check:
            try:
                atom1 = int(atom1)
                atom2 = int(atom2)
            except Exception as e:
                raise TypeError("atom numbers cannot be cast to int!")

            self._check_atom_number(atom1)
            self._check_atom_number(atom2)

        return _dist(self.get_vector(atom1, check=False), self.get_vector(atom2, check=False))

    def get_sq_distance(self, atom1, atom2, check=True):
        """
        Wrapper to compute squared distance between two atoms -- optimized for speed!

        Args:
            atom1 (int): number of the first atom
            atom2 (int): number of the second atom
            check (Bool): whether to validate input data (can be overridden to prevent slow double-checking)

        Returns:
            the squared distance
        """
        if check:
            try:
                atom1 = int(atom1)
                atom2 = int(atom2)
            except Exception as e:
                raise TypeError("atom numbers cannot be cast to int!")

            self._check_atom_number(atom1)
            self._check_atom_number(atom2)

        return np.sum(np.square(self.get_vector(atom1, atom2, check=False)))

    def get_angle(self, atom1=None, atom2=None, atom3=None, check=True, _angle=compute_angle_between, atoms=None):
        """
        Wrapper to compute angle between three atoms.

        This function is relatively slow (rate-limiting for certain applications), so performance boosts have been implemented (e.g. preloading ``_angle``).

        Args:
            atom1 (int): number of the first atom
            atom2 (int): number of the second atom
            atom3 (int): number of the third atom
            check (Bool): whether to validate input data (can be overridden to prevent slow double-checking)
            _angle (function): function usd to compute angle
            atoms (list): list of atom numbers

        Returns:
            the angle, in degrees
        """
        if (atom1 is None) and (atom2 is None) and (atom3 is None):
            assert isinstance(atoms, (list, np.ndarray)), "atom numbers need to come from fields or list!"
            assert len(atoms) == 3, "need 3 atom numbers to get angle"
            atom1 = atoms[0]
            atom2 = atoms[1]
            atom3 = atoms[2]

        if check:
            try:
                atom1 = int(atom1)
                atom2 = int(atom2)
                atom3 = int(atom3)
            except Exception as e:
                raise TypeError("atom numbers cannot be cast to int!")

            self._check_atom_number(atom1)
            self._check_atom_number(atom2)
            self._check_atom_number(atom3)

        v1 = self.get_vector(atom1, check=False)
        v2 = self.get_vector(atom2, check=False)
        v3 = self.get_vector(atom3, check=False)

        return _angle(v1 - v2, v3 - v2)

    def get_dihedral(self, atom1=None, atom2=None, atom3=None, atom4=None, check=True, _dihedral=compute_dihedral_between, atoms=None):
        """
        Wrapper to compute dihedral angle between four atoms.

        This function is relatively slow (rate-limiting for certain applications), so performance boosts have been implemented (e.g. preloading ``_dihedral``).

        Args:
            atom1 (int): number of the first atom
            atom2 (int): number of the second atom
            atom3 (int): number of the third atom
            atom4 (int): number of the fourth atom
            check (Bool): whether to validate input data (can be overridden to prevent slow double-checking)
            _dihedral (function): function used to compute dihedral
            atoms (list): list of atom numbers

        Returns:
            the dihedral angle, in degrees
        """
        if (atom1 is None) and (atom2 is None) and (atom3 is None) and (atom4 is None):
            assert isinstance(atoms, (list, np.ndarray)), "atom numbers need to come from fields or list!"
            assert len(atoms) == 4, "need 4 atom numbers to get dihedral angle"
            atom1 = atoms[0]
            atom2 = atoms[1]
            atom3 = atoms[2]
            atom4 = atoms[3]

        if check:
            try:
                atom1 = int(atom1)
                atom2 = int(atom2)
                atom3 = int(atom3)
                atom4 = int(atom4)
            except Exception as e:
                raise TypeError("atom numbers cannot be cast to int!")

            self._check_atom_number(atom1)
            self._check_atom_number(atom2)
            self._check_atom_number(atom3)
            self._check_atom_number(atom4)

        return _dihedral(
            self.get_vector(atom1, check=False),
            self.get_vector(atom2, check=False),
            self.get_vector(atom3, check=False),
            self.get_vector(atom4, check=False),
        )

    def get_bond_order(self, atom1, atom2):
        """
        Wrapper to get bond order between two atoms.

        Args:
            atom1 (int): number of the first atom
            atom2 (int): number of the second atom

        Returns:
            the bond order
        """
        self._check_atom_number(atom1)
        self._check_atom_number(atom2)

        if self.bonds.has_edge(atom1, atom2):
            return self.bonds[atom1][atom2]["weight"]
        else:
            return 0

    def are_connected(self, atom1, atom2):
        """
        Wrapper to tell if two atoms are connected.
        """
        self._check_atom_number(atom1)
        self._check_atom_number(atom2)

        if atom1 in self._get_fragment_containing(atom2):
            return True
        else:
            return False

    def get_atoms_by_symbol(self, symbol):
        """
        Returns all the numbers of atoms of type ``symbol`` in the molecule.
        """
        if not isinstance(symbol, str):
            raise TypeError("symbol {symbol} must be a string")

        number = get_number(symbol)
        atoms = []

        for index, atom in enumerate(self.atomic_numbers, start=1):
            if atom == number:
                atoms.append(index)

        return atoms

    def get_heavy_atoms(self):
        """
        Returns a list of all the heavy atoms in the molecule (i.e., not hydrogen), for array indexing.
        """
        atoms = []

        for index, atom in enumerate(self.atomic_numbers, start=1):
            if atom != 1:
                atoms.append(index)

        return atoms

    def get_adjacent_atoms(self, atom):
        """
        Returns a list of the neighbors of ``atom``. If ``atom`` has no neighbors, an empty list will be returned.
        """
        try:
            atom = int(atom)
        except Exception as e:
            raise TypeError(f"atom number {atom} cannot be cast to int!")

        self._check_atom_number(atom)

        return list(self.bonds.neighbors(atom))

    def num_atoms(self):
        return len(self.atomic_numbers)

    def rms_distance_between_atoms(self):
        """
        Returns the RMS distance (in Angstroms) between every pair of atoms - a quick, easy-to-calculate proxy for minimizing steric clashes.
        """
        distance = 0
        for i in range(1, self.num_atoms() + 1):
            for j in range(1, self.num_atoms() + 1):
                if i == j:
                    continue
                distance += self.get_distance(i, j) ** 2

        return math.sqrt(distance) / self.num_atoms()

    def optimize_dihedral(self, atom1, atom2, atom3, atom4, step=10):
        """
        Minimizes the value of ``self.rms_distance_between_atoms`` for the given dihedral, in one-degree increments.
        A cheap alternative to geometry optimization using *ab initio* methods or density functional theory.

        Args:
            atom1 (int): atom number of the first atom in the dihedral
            atom2 (int): atom number of the second atom in the dihedral
            atom3 (int): atom number of the third atom in the dihedral
            atom4 (int): atom number of the fourth atom in the dihedral
            step (float): explore angles from 0 to 360 with this stepsize in degrees

        Returns:
            the final value of the angle
        """
        self._check_atom_number(atom1)
        self._check_atom_number(atom2)
        self._check_atom_number(atom3)
        self._check_atom_number(atom4)

        best_angle = 0
        best_dist = 0

        for angle in range(0, 360, step):
            self.set_dihedral(atom1, atom2, atom3, atom4, angle)
            if self.rms_distance_between_atoms() > best_dist:
                best_dist = self.rms_distance_between_atoms()
                best_angle = angle

        self.set_dihedral(atom1, atom2, atom3, atom4, best_angle)
        return best_angle

    def atom_string(self, atom):
        """
        Returns the elemental symbol and the atom number for a given atom.

        For example, ``methane.atom_string(1)`` might return "C1".

        Args:
            atom (int): number of the atom

        Returns:
            the aforementioned atom string
        """
        try:
            atom = int(atom)
        except Exception as e:
            raise ValueError("atom cannot be cast to int")

        self._check_atom_number(atom)

        return f"{get_symbol(self.atomic_numbers[atom])}{atom}"

    def perturb(self, size=0.005):
        """
        This function can be used to generate a slightly different molecule in cases where numerical (or geometric) converge is problematic.

        It adds a random variable (sampled from a normal distribution, centered at 0 with stddev ``size`) to every number in ``self.geometry``.

        Args:
            size (float): stddev of the normal distribution

        Returns:
            the Molecule object
        """
        geometry = self.geometry
        random = np.random.normal(scale=size, size=geometry.shape)

        self.geometry = geometry + random
        return self

    def center(self):
        """
        Moves the centroid to the origin.
        """
        atoms = np.arange(1, self.num_atoms()+1)
        self.translate_molecule(-self.geometry[atoms].mean(axis=0))
        return self

    @classmethod
    def combine_molecules(cls, molecule1, molecule2):
        """
        Combine two molecules into one final molecule.

        Bonding information is not currently preserved.

        Args:
            molecule1 (Molecule): 1st molecule
            molecule2 (Molecule): 2nd molecule

        Returns:
            new ``Molecule`` object
        """

        atoms = np.hstack((molecule1.atomic_numbers.T, molecule2.atomic_numbers.T)).view(cctk.OneIndexedArray)
        geoms = np.vstack((molecule1.geometry, molecule2.geometry)).view(cctk.OneIndexedArray)
        charge = molecule1.charge + molecule2.charge

        s1 = (molecule1.multiplicity - 1) / 2
        s2 = (molecule2.multiplicity - 1) / 2
        multiplicity = (s1+s2) * 2 + 1

        return Molecule(atoms, geoms, charge=charge, multiplicity=multiplicity)

    def volume(self, pts_per_angstrom=10, qhull=False):
        """
        Returns volume calculated using the Gavezotti algorithm (JACS, 1983, 105, 5220). Relatively slow.
        If MemoryError, defaults to a qhull-based approach (accurate in the limit as number of atoms goes to infinity)

        Args:
            pts_per_angstrom (int): how many grid points to use per Å - time scales as O(n**3) so be careful!
            qhull (bool): use faster QHull algorithm

        Returns:
            volume in Å**3
        """
        if not qhull:
            try:
                assert isinstance(pts_per_angstrom, int), "Need an integer number of pts per Å!"
                assert pts_per_angstrom > 0, "Need a positive integer of pts per Å!"

                box_max = np.max(self.geometry.view(np.ndarray), axis=0) + 4
                box_min = np.min(self.geometry.view(np.ndarray), axis=0) - 4

                box_volume = (box_max[0] - box_min[0]) * (box_max[1] - box_min[1]) * (box_max[2] - box_min[2])

                x_vals = np.linspace(box_min[0], box_max[0], int((box_max[0] - box_min[0]) * pts_per_angstrom))
                y_vals = np.linspace(box_min[1], box_max[1], int((box_max[1] - box_min[1]) * pts_per_angstrom))
                z_vals = np.linspace(box_min[2], box_max[2], int((box_max[2] - box_min[2]) * pts_per_angstrom))

                # h4ck3r
                box_pts = np.stack([np.ravel(a) for a in np.meshgrid(x_vals, y_vals, z_vals)], axis=-1)

                # caching to speed call
                vdw_radii = {z: get_vdw_radius(z) for z in set(self.atomic_numbers)}
                radii_per_atom = np.array([vdw_radii[z] for z in self.atomic_numbers]).reshape(-1,1)

                # this is the slow part since it's approximately a zillion operations
                dists_per_atom = cdist(self.geometry.view(np.ndarray), box_pts)
                occupied = np.sum(np.max(dists_per_atom < radii_per_atom, axis=0))

                percent_occupied = occupied / box_pts.shape[0]
                return percent_occupied * box_volume
            except MemoryError:
                qhull = True

        if qhull:
            import scipy
            hull = scipy.spatial.ConvexHull(self.geometry.view(np.ndarray))
            return hull.volume

    def swap_atom_numbers(self, atom1, atom2):
        """
        Interchanges the numbers of ``atom1`` and ``atom2``.

        Args:
            atom1 (int): number of 1st atom
            atom2 (int): number of 2nd atom

        Returns
            new ``Molecule`` object (does not modify in-place)
        """
        self._check_atom_number(atom1)
        self._check_atom_number(atom2)
        mol = copy.deepcopy(self)

        z1 = mol.atomic_numbers[atom1]
        z2 = mol.atomic_numbers[atom2]
        g1 = copy.deepcopy(mol.geometry[atom1])
        g2 = copy.deepcopy(mol.geometry[atom2])

        mol.atomic_numbers[atom2] = z1
        mol.atomic_numbers[atom1] = z2
        mol.geometry[atom2] = g1
        mol.geometry[atom1] = g2

        mapping = {atom2: atom1, atom1: atom2}
        mol.bonds = nx.relabel_nodes(mol.bonds, mapping, copy=True)
        return mol

    def epimerize(self, center_atom, substituent1, substituent2):
        """
        Epimerizes ``center_atom`` by exchanging the groups corresponding to ``substituent1`` and ``substituent2``.
        Both substituents must be bonded to the center atom!

        Args:
            center_atom (int): number of middle atom
            substituent1 (int): number of 1st atom
            substituent1 (int): number of 2nd atom

        Returns
            new ``Molecule`` object (does not modify in-place)
        """

        self._check_atom_number(center_atom)
        self._check_atom_number(substituent1)
        self._check_atom_number(substituent2)

        assert self.bonds.number_of_edges() > 0, "need a bond graph to perform this operation -- try calling self.assign_connectivity()!"

        adj = self.get_adjacent_atoms(center_atom)
        assert len(adj) == 4, "center atom must be making 4 bonds!"
        assert substituent1 in adj, "1st substituent is not bonded to center atom!"
        assert substituent2 in adj, "2nd substituent is not bonded to center atom!"

        #### remove both substituents
        mol, group1, mmap1, gmap1  = cctk.Group.remove_group_from_molecule(self, center_atom, substituent1, return_mapping=True)
        mol, group2, mmap2, gmap2  = cctk.Group.remove_group_from_molecule(mol, mmap1[center_atom], mmap1[substituent2], return_mapping=True)

        h1 = mol.num_atoms() - 1
        h2 = mol.num_atoms()

        #### add them back in the opposite fashion
        mol, mmap3, gmap3 =  cctk.Group.add_group_to_molecule(mol, group2, h1, return_mapping=True)
        mol = cctk.Group.add_group_to_molecule(mol, group1, mmap3[h2])

        #### relabel new graph to match original molecule
        which = top.get_stereogenic_centers(self)
        which.remove(center_atom)
        return mol.renumber_to_match(self, check_chirality=which)

    def renumber_to_match(self, model, check_chirality="all"):
        """
        Renumbers atoms to match ``model`` (must have isomorphic bond graph). Returns a copy of ``self`` with renumbered atoms.

        Args:
            model (cctk.Molecule): isomorphic molecule to renumber by
            check_chirality (list of atomic numbers): atomic numbers to check, to prevent inversion due to graph isomorphism.
                Alternatively ``None`` will prevent any checking and "all" will use ``cctk.topology.get_exchangable_centers()``.

        Returns:
            new ``Molecule`` object
        """

        assert self.bonds.number_of_edges() > 0, "need a bond graph to perform this operation -- try calling self.assign_connectivity()!"

        #### use networkx to generate mapping
        #### you need the node matcher to distinguish between e.g. H, F, Cl
        self._add_atomic_numbers_to_nodes()
        model._add_atomic_numbers_to_nodes()
        nm = nx.algorithms.isomorphism.categorical_node_match("atomic_number", 0)

        match = nx.algorithms.isomorphism.GraphMatcher(model.bonds, self.bonds, node_match=nm)
        assert match.is_isomorphic(), "can't renumber non-isomorphic graphs!"
        new_ordering = [match.mapping[x] for x in range(1, self.num_atoms() + 1)]
        inv_mapping = {v:k  for k,v in match.mapping.items()} # bit kludgy but works

        #### create renumbered molecule
        mol = copy.deepcopy(self)
        mol.atomic_numbers = self.atomic_numbers[new_ordering]
        mol.geometry = self.geometry[new_ordering]
        mol.bonds = nx.relabel_nodes(self.bonds, mapping=inv_mapping, copy=True)

        if check_chirality == "all":
            check_chirality = top.get_exchangeable_centers(mol)

        #### diastereotopic protons get scrambled by the above code so we gotta go through and fix all of them
        #### this happens because networkx doesn't store chirality - a known limitation of graph molecular encoding!
        if isinstance(check_chirality, list):
            #### find all the differences and exchange them
            model_report = top.get_chirality_report(model, check_chirality)

            #### generate all meso ring permutations
            candidates = top.flip_meso_rings(mol, atoms=check_chirality)

            #### for each, try flipping configuration of all centers
            for candidate in candidates:
                report = top.get_chirality_report(candidate, check_chirality)
                for center in check_chirality:
                    if model_report[center] != report[center]:
                        try:
                            candidate = top.exchange_identical_substituents(candidate, center)
                        except ValueError as e:
                            break

                #### check that we actually fixed all the problems
                mol_report = top.get_chirality_report(candidate, check_chirality)
                all_good = True
                for center in check_chirality:
                    if mol_report[center] != model_report[center]:
                        all_good = False
                        break
                #### if we did, then return
                if all_good:
                    return candidate

        raise ValueError("can't get a proper renumbering: are you *sure* these two molecules can have the same chirality?")

    def _add_atomic_numbers_to_nodes(self):
        """
        Add the atomic numbers to each node attribute, to allow for distinguishment of F and H during graph renumbering.
        """
        nx.set_node_attributes(self.bonds, {z: {"atomic_number": self.atomic_numbers[z]} for z in range(1, self.num_atoms() +  1)})

    def is_atom_in_ring(self, atom):
        assert self.bonds.number_of_edges() > 0, "need a bond graph to perform this operation -- try calling self.assign_connectivity()!"
        cycles = nx.cycle_basis(self.bonds, root=atom)
        for cycle in cycles:
            if atom in cycle:
                return True
        return False

    def get_components(self):
        """
        Returns a list of all the connected components in a molecule.
        """
        assert self.bonds.number_of_edges() > 0, "need a bond graph to perform this operation -- try calling self.assign_connectivity()!"
        fragments = nx.connected_components(self.bonds)
        return [list(f) for f in list(fragments)]

    def limit_solvent_shell(self, solute=0, num_atoms=0, num_solvents=10, distance_from_atom=None, return_idxs=False):
        """
        Automatically detects solvent molecules and removes them until you have a set number of solvents or atoms.

        The "distance" between molecules is the minimum of the pairwise atomic distances, to emphasize inner-sphere interactions.

        Args:
            solute (int): which fragment is the solute, 0-indexed
            num_atoms (int): remove atoms until there are this number (modulo the size of a solvent molecule)
            num_solvents (int): remove solvent molecules until there are this number
            distance_from_atom (int): if you want to find molecules closest to a given atom in the solute, specify the atom number here.
                if this atom is not in the solute fragment, an exception will be raised.
            return_idxs (bool): if True, indices of atoms that would be in the new molecule are returned. no change is made to ``self``.

        Returns:
            new ``Molecule`` object
        """
        assert isinstance(num_atoms, int)
        assert isinstance(num_solvents, int)

        fragments = self.get_components()
        solute_x = self.geometry[fragments[solute]].view(np.ndarray)

        if distance_from_atom:
            assert distance_from_atom in fragments[solute], f"{distance_from_atom} is not in the solute fragment"
            solute_x = self.geometry[[distance_from_atom]].view(np.ndarray)

        distances = np.zeros(shape=len(fragments))
        for i, f in enumerate(fragments):
            if i == solute:
                distances[i] = 0
            else:
                solvent_x = self.geometry[f].view(np.ndarray)
                # cdist is absurdly fast
                pairwise_distances = cdist(solvent_x, solute_x)
                distances[i] = np.min(pairwise_distances)

        mol = copy.deepcopy(self)

        #### reverse order - farthest away comes first
        order = np.argsort(distances)[::-1]

        current_num_solvents = len(fragments) - 1
        current_num_atoms = mol.num_atoms()

        to_remove = []
        for i in order:
            for j in fragments[i]:
                to_remove.append(j)
                current_num_atoms += -1
            current_num_solvents += -1

            if current_num_atoms <= num_atoms or num_solvents == current_num_solvents:
                if return_idxs:
                    all_idxs = set(range(1,self.num_atoms()))
                    return list(all_idxs - set(to_remove))
                else:
                    #### have to remove in reverse direction for indexing consistency
                    for j in sorted(to_remove, reverse=True):
                        mol.remove_atom(j)
                    return mol

    def center_periodic(self, center, side_length):
        """
        Adjusts a molecule to be in the center of a cube, moving all other molecules accordingly. Bonded subgroups will be moved as a unit.

        For analysis of MD files with periodic boundary conditions.

        Args:
            center (int): atomic number to center
            side_length (float): length of side, in Å
        """
        self._check_atom_number(center)
        assert isinstance(side_length, (int, float))
        assert side_length > 0

        #### Center the atom of interest
        self.geometry += -1 * self.geometry[center]
        self.geometry += side_length / 2

        for f in self.get_components():
            centroid = np.mean(self.geometry[f], axis=0)
            self.geometry[f] += -1 * np.floor_divide(centroid, side_length) * side_length

        return self

    @classmethod
    def new_from_name(cls, name):
        """
        Create a new ``Molecule`` instance using ``rdkit``.
        """
        assert isinstance(name, str)
        from urllib.request import urlopen

        try:
            url_name = re.sub(" ", "%20", name)
            url = 'http://cactus.nci.nih.gov/chemical/structure/' + url_name + '/smiles'
            smiles = urlopen(url, timeout=5).read().decode('utf8')
            return cls.new_from_smiles(smiles)
        except Exception as e:
            raise ValueError(f"something went wrong auto-generating molecule {name}:\nurl: {url}\n{e}")

    @classmethod
    def new_from_smiles(cls, smiles):
        """
        Create a new ``Molecule`` instance using ``rdkit``.
        """
        assert isinstance(smiles, str)

        try:
            from rdkit.Chem import AllChem as Chem
        except ImportError as e:
            raise ImportError(f"``rdkit`` must be installed for this function to work!\n{e}")

        try:
            rdkm = Chem.MolFromSmiles(smiles)
            rdkm = Chem.AddHs(rdkm)
            Chem.EmbedMolecule(rdkm)
            Chem.MMFFOptimizeMolecule(rdkm)

            nums = []
            for atom in rdkm.GetAtoms():
                nums.append(atom.GetAtomicNum())
            geom = rdkm.GetConformers()[0].GetPositions()

            return cls(nums, geom)

        except Exception as e:
            raise ValueError(f"something went wrong auto-generating molecule {smiles}:\n{e}")

    def fragment(self):
        """
        Returns list of ``cctk.Molecule`` objects based on the bond-connected components of ``self``.
        """
        fragments = list()
        indices = self.get_components()
        for idx in indices:
            mol = cctk.Molecule(self.atomic_numbers[idx], self.geometry[idx]).assign_connectivity()
            fragments.append(mol)
        return fragments

    def get_symmetric_atoms(self):
        """
        Returns lists of symmetric atoms, as defined in ``cctk.load_group``.

        Useful for NMR spectroscopy, etc.
        """
        from cctk.load_groups import group_iterator

        symmetric_sets = []
        for group in group_iterator(symmetric_only=True):
            # this gives us a list of dictionaries mapping from self.atomic_numbers to group numbers
            matches = top.find_group(self, group)

            for m in matches:
                i = {v: k for k,v in m.items()}
                for n in group.isomorphic:
                    symmetric_sets.append([i[idx] for idx in n])

        #### some groups overlap (e.g. methyl and t-butyl), so now we collapse the overlapping sets
        for i, s1 in enumerate(symmetric_sets):
            for j, s2 in enumerate(symmetric_sets[i+1:]):
                if set(s1).intersection(set(s2)):
                    symmetric_sets[i + j + 1] = list(set(s1).union(s2))
                    symmetric_sets[i] = None # can't delete yet - messes up indexing

        #### now we delete
        symmetric_sets = list(filter(None, symmetric_sets))
        return symmetric_sets

    def atomic_symbols(self):
        """
        Return list of atomic symbols.
        """
        symbols = {z: get_symbol(z) for z in set(self.atomic_numbers)}
        return [symbols[z] for z in self.atomic_numbers]

    def optimize(self, inplace=True, nprocs=1, return_energy=False):
        """
        Optimize molecule at the GFN2-xtb level of theory.

        Args:
            inplace (Bool): whether or not to return a new molecule or simply modify ``self.geometry``
            nprocs (int): number of processors to use
            return_energy (Bool): whether to return energy or not
        """
        import cctk.optimize as opt
        assert isinstance(nprocs, int), "nprocs must be int!"
        optimized, energy = opt.optimize_molecule(self, nprocs=nprocs, return_energy=True)

        if inplace:
            self.geometry = optimized.geometry
            if return_energy:
                return self, energy
            else:
                return self
        else:
            if return_energy:
                return optimized, energy
            else:
                return optimized

    def compute_energy(self, nprocs=1):
        """
        Compute energy of molecule at the GFN2-xtb level of theory.

        Args:
            nprocs (int): number of processors to use
        """
        import cctk.optimize as opt
        assert isinstance(nprocs, int), "nprocs must be int!"
        energy = opt.get_energy(self, nprocs=nprocs)
        return energy

    def csearch(self, nprocs=1, constraints=[], logfile=None, noncovalent=False, use_tempdir=True, gfn=2, additional_flags=None):
        """
        Optimize molecule at the GFN2-xtb level of theory.

        Args:
            nprocs (int): number of processors to use
            constraints (list): atoms numbers to freeze
            noncovalent (bool): whether or not to use non-covalent settings
            logfile (str): file to write ongoing ``crest`` output to
            use_tempdir (bool): write intermediate files to hidden directory (as opposed to current directory)
            gfn (int or ``ff``): level of theory, either 1, 2, or ``ff``
            additional_flags (str): additional flags for command line

        Returns
            ConformationalEnsemble
        """
        import cctk.optimize as opt
        assert isinstance(nprocs, int), "nprocs must be int!"
        return opt.csearch(molecule=self, nprocs=nprocs, constraints=constraints, noncovalent=noncovalent, logfile=logfile, use_tempdir=use_tempdir, gfn=gfn, additional_flags=additional_flags)

    def num_neighbors_by_atom(self):
        """
        Returns a list of the number of neighbors of each atom.
        """
        result = []
        for i in range(self.num_atoms()):
            result.append(len(self.get_adjacent_atoms(i)))
        return result

    def atoms_moving_in_imaginary(self, max_num=5, percent_cutoff=0.03, return_string=False):
        """
        Returns atoms moving in imaginary, ranked by how much they're moving.

        Args:
            max_num (int): how many atoms max to return
            percent_cutoff (float): threshold for what percent of total TS movement qualifies as "movement"
            return_string (bool): whether or not to return a formatted string report

        Returns:
            list of atomic numbers or string
        """
        imaginary = 0
        ts_mode = None
        for mode in self.vibrational_modes:
            if mode.frequency < imaginary:
                imaginary = mode.frequency
                ts_mode = mode

        if ts_mode is None:
            if return_string:
                return ""
            else:
                return None

        displacements = np.linalg.norm(ts_mode.displacements.view(np.ndarray), axis=-1)

        atoms_ranked = np.argsort(displacements)[::-1] + 1
        percent_movement = np.sort(displacements)[::-1] / np.sum(displacements)

        return_list, string = list(), ""
        for atom, percent in zip(atoms_ranked, percent_movement):
            if percent > percent_cutoff and len(return_list) <= max_num:
                return_list.append(atom)
                string += f"{self.atom_string(atom)} ({percent:.1%}), "
            else:
                if return_string:
                    return string[:-2]
                else:
                    return return_list


    def to_string(self):
        """
        Save the current molecule as a string, for subsequent loading. Not human-readable.

        Vibrational modes are currently not saved.
        """
        # name, charge, multiplicity need no encoding
        atomic_number_encoding = numpy_to_bytes(self.atomic_numbers.view(np.ndarray))
        geometry_encoding = numpy_to_bytes(self.geometry.view(np.ndarray))
        bonds_encoding = numpy_to_bytes(nx.convert_matrix.to_numpy_array(self.bonds))

        if self.name is None:
            self.name = "name"

        cctk_version = pkg_resources.get_distribution("cctk").version

        store_dict = {
            "name": self.name,
            "charge": self.charge,
            "multiplicity": self.multiplicity,
            "atomic_numbers": atomic_number_encoding,
            "geometry": geometry_encoding,
            "bonds": bonds_encoding,
            "cctk_version": cctk_version,
        }

        return yaml.dump(store_dict)

    @classmethod
    def from_string(cls, string, check_version=True):
        """
        Loads a ``cctk.Molecule`` object from a string.

        Arguments:
            string (str): stringified version of the molecule
            check_version (bool): whether version consistency should be enforced
        """

        try:
            store_dict = yaml.safe_load(string)

            if check_version:
                cctk_version = pkg_resources.get_distribution("cctk").version
                assert cctk_version == store_dict["cctk_version"], f"Warning: the data was saved in cctk {store_dict['cctk_version']} but is being loaded in cctk {cctk_version}!"

            atomic_numbers = bytes_to_numpy(store_dict["atomic_numbers"]).astype(np.int8)
            geometry = bytes_to_numpy(store_dict["geometry"]).astype(np.float32)
            bonds = nx.convert_matrix.from_numpy_array(bytes_to_numpy(store_dict["bonds"]))

            mol = cls(
                atomic_numbers,
                geometry,
                bonds=bonds,
                charge=store_dict["charge"],
                multiplicity=store_dict["multiplicity"],
                name=store_dict["name"],
                checks=False, # trust nx data implicitly
            )

            return mol

        except Exception as e:
            raise ValueError(f"this stringified Molecule fails import: {e}")

    def coulomb_analysis(self, atoms1, atoms2, charges):
        """
        Computes the net Coulomb forces between atoms ``atoms1`` and atoms ``atoms2``.
        """
        if isinstance(charges, np.ndarray):
            charges = charges.view(cctk.OneIndexedArray)
        elif isinstance(charges, list):
            charges = cctk.OneIndexedArray(charges)

        assert isinstance(charges, cctk.OneIndexedArray), "charges must be cctk.OneIndexedArray"
        assert len(charges) == self.num_atoms(), "need a charge for every atom"
        assert isinstance(atoms1, list)
        assert isinstance(atoms2, list)

        q1 = charges[atoms1]
        q2 = charges[atoms2]

        # need to convert to Bohr
        r1 = self.geometry[atoms1] / 0.529
        r2 = self.geometry[atoms2] / 0.529

        R = cdist(r1, r2)**2
        Q = np.outer(q1, q2)

        energy = 0
        for i in range(len(atoms1)):
            assert atoms1[i] not in atoms2, "lists must be non-overlapping"
            for j in range(len(atoms2)):
                energy += Q[i][j] / R[i][j]

        return energy * 627.509 # convert to kcal/mol
