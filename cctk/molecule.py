import sys, re, math
import numpy as np
import networkx as nx

from functools import lru_cache

from cctk import OneIndexedArray
from cctk.helper_functions import (
    get_symbol,
    get_number,
    compute_rotation_matrix,
    compute_distance_between,
    compute_angle_between,
    compute_dihedral_between,
    compute_unit_vector,
    get_covalent_radius,
    get_isotopic_distribution,
)

class Molecule:
    """
    Class that represents a single molecule, abstractly.

    In contrast to typical Python behavior, ``atomic_numbers`` and ``geometry`` are indexed from one, to simplify interfacing with computational chemistry programs.
    This has been done by defining a custom wrapper for ``numpy.ndarray`` called ``cctk.OneIndexedArray``.

    All other datatypes are indexed from 0.

    Attributes:
        name (str): for identification, optional
        atomic_numbers (cctk.OneIndexedArray, dtype=np.int8): list of atomic numbers
        geometry (cctk.OneIndexedArray): list of 3-tuples of xyz coordinates - same ordering as ``atomic_numbers``
        bonds (nx.Graph): Graph object containing connectivity information (1-indexed)
        charge (int): the charge of the molecule
        multiplicity (int): the spin state of the molecule (1 corresponds to singlet, 2 to doublet, 3 to triplet, etc. -- so a multiplicity of 1 is equivalent to S=0)
    """

    def __init__(self, atomic_numbers, geometry, name=None, bonds=None, charge=0, multiplicity=1):
        """
        Create new Molecule object, and assign connectivity if needed.

        ``bonds`` must be a list of edges (i.e. an n x 2 ``numpy`` array).
        """
        if len(atomic_numbers) != len(geometry):
            raise ValueError("length of geometry and atomic_numbers does not match!")

        if not all(isinstance(z, np.int8) for z in atomic_numbers) or atomic_numbers.size == 0:
            raise ValueError("invalid atom list")

        if len(geometry) == 0:
            raise ValueError("invalid geometry list")

        try:
            geometry = np.array(geometry)
            geometry = geometry.astype(float)
        except:
            raise TypeError("geometry cannot be cast to ndarray of floats!")

        if not all(all(isinstance(y, float) for y in x) for x in geometry):
            raise TypeError("each element of self.geometry must be a 3-tuple")

        if name and not isinstance(name, str):
            raise TypeError("name must be a string!")

        for atom in atomic_numbers:
            try:
                get_symbol(atom)
            except ValueError:
                raise ValueError(f"unknwon atomic number {atom}")

        if not isinstance(charge, int):
            try:
                charge = int(charge)
            except:
                raise TypeError("charge must be integer or castable to integer!")

        if not isinstance(multiplicity, int):
            try:
                multiplicity = int(multiplicity)
            except:
                raise TypeError("multiplicity must be positive integer or castable to positive integer")
        assert multiplicity > 0, "multiplicity must be positive"

        self.atomic_numbers = OneIndexedArray(atomic_numbers, dtype=np.int8)
        self.geometry = OneIndexedArray(geometry)

        if bonds:
            for bond in bonds:
                if len(bond) != 2:
                    raise ValueError("while 3-center bonding is possible, it's a no-go in cctk")
                self._check_atom_number(bond[0])
                self._check_atom_number(bond[1])

        self.name = name
        self.multiplicity = multiplicity
        self.charge = charge

        self.bonds = nx.Graph()
        self.bonds.add_nodes_from(range(1, len(atomic_numbers) + 1))
        if bonds:
            for bond in bonds:
                self.add_bond(bond[0], bond[1])

    def assign_connectivity(self, cutoff=0.5):
        """
        Automatically recalculates bonds based on covalent radii. If two atoms are closer than the sum of their covalent radii + 0.5 Angstroms, then they are considered bonded.

        Args:
            cutoff (float): the threshold (in Angstroms) for how close two covalent radii must be to be considered bonded

        Returns:
            self
        """

        for i in range(1, self.num_atoms() + 1):
            for j in range(i + 1, self.num_atoms() + 1):
                distance = self.get_distance(i, j)
                r_i = get_covalent_radius(self.get_atomic_number(i))
                r_j = get_covalent_radius(self.get_atomic_number(j))

                # 0.5 A distance is used by RasMol and Chime (documentation available online) and works well, empirically
                if distance < (r_i + r_j + cutoff):
                    self.add_bond(i, j)
                elif self.get_bond_order(i, j):
                    self.remove_bond(i, j)

        return self

    def check_for_conflicts(self, min_buffer=-1, group1=None, group2=None):
        """
        Automatically checks for conflicts based on covalent radii. If two atoms are closer than the sum of their covalent radii + buffer, then they are considered clashing.
        If `group1` and `group2` are selected, then conflicts will only be evaluated between these two groups of atoms.

        Args:
            min_buffer (float): the threshold (in Angstroms) for how close two covalent radii must be to be considered clashing. -1.0 A is default, for no particular reason.
            group1 (list): atoms to evaluate against `group2` (if `None`, defaults to all atoms)
            group2 (list): atoms to evaluate against `group1` (if `None`, defaults to all atoms)

        Returns:
            True if there are no conflicts, False (+ error) if there are
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
                if distance < (r_i + r_j + min_buffer):
                    raise ValueError(f"atoms {i} and {j} are too close - distance {distance} A!")
                    return False

        return True

    def add_bond(self, atom1, atom2, bond_order=1):
        """
        Adds a new bond to the bond graph, or updates the existing bond order. Will not throw an error if the bond already exists.

        Args:
            atom1 (int): the number of the first atom
            atom2 (int): the number of the second atom
            bond_order (int): bond order of bond between atom1 and atom2
        """
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
        if not isinstance(number, int):
            raise TypeError("atom number must be integer")

        if number > self.num_atoms():
            raise ValueError(f"atom number {number} too large!")

        if number <= 0:
            raise ValueError(f"atom number {number} invalid: must be a positive integer!")

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
            if "H" in elements:
                elements.remove("H")
                formula += f"H{formula_dict['H']}"

            if "C" in elements:
                elements.remove("C")
                formula += f"C{formula_dict['C']}"

            for element in sorted(elements):
                formula += f"{element}{formula_dict[element]}"

            return formula

    #### very fast but causes errors sometimes... so i'm commenting this out until further consultation.
    #    @lru_cache(maxsize=32)
    def _get_bond_fragments(self, atom1, atom2, bond_order=1):
        """
        Returns the pieces of a molecule that one would obtain by breaking the bond between two atoms. Will throw ``ValueError`` if the atoms are in a ring.
        Useful for distance/angle/dihedral scans -- one fragment can be moved and the other held constant.

        Args:
            atom1 (int): the number of the first atom
            atom2 (int): the number of the second atom
            bond_order (int): bond order of bond between atom1 and atom2

        Returns:
            fragment1: the list of atoms in fragment 1 (containing atom1)
            fragment2: the list of atoms in fragment 2 (containing atom2)

        """

        self._check_atom_number(atom1)
        self._check_atom_number(atom2)

        if (not isinstance(bond_order, int)) or (bond_order < 0):
            raise ValueError("invalid bond order!")

        if self.bonds.has_edge(atom1, atom2):
            self.bonds.remove_edge(atom1, atom2)

            fragments = nx.connected_components(self.bonds)
            fragment1 = []
            fragment2 = []

            for fragment in fragments:
                if atom1 in fragment:
                    if atom2 in fragment:
                        raise ValueError(f"Atom {atom1} and atom {atom2} are in a ring or otherwise connected!")
                    else:
                        fragment1 = fragment
                if atom2 in fragment:
                    fragment2 = fragment

            self.bonds.add_edge(atom1, atom2, weight=bond_order)
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
                break

    def set_distance(self, atom1, atom2, distance, move="group"):
        """
        Adjusts the ``atom1`` -- ``atom2`` bond length to be a fixed distance by moving atom2.

        If ``move`` is set to "group", then all atoms bonded to ``atom2`` will also be moved.

        If ``move`` is set to "atom", then only atom2 will be moved.

        Args:
            atom1 (int): the number of the first atom
            atom2 (int): the number of the second atom
            distance (float): distance in Angstroms of the final bond
            move (str): determines how fragment moving is handled

        Returns:
            the Molecule object
        """

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

    def set_angle(self, atom1, atom2, atom3, angle, move="group"):
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

        Returns:
            the Molecule object
        """

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
        except:
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

    def set_dihedral(self, atom1, atom2, atom3, atom4, dihedral, move="group34", check_result=True):
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

        Returns:
            the Molecule object
        """

        self._check_atom_number(atom1)
        self._check_atom_number(atom2)
        self._check_atom_number(atom3)
        self._check_atom_number(atom4)

        for x in [atom1, atom2, atom3, atom4]:
            for y in [atom1, atom2, atom3, atom4]:
                if x <= y:
                    continue
                else:
                    if self.get_sq_distance(x, y, check=False) < 0.001:
                        raise ValueError(f"atom {x} and atom {y} are too close!")

        try:
            dihedral = float(dihedral)
        except:
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
            return

        #### now the real work begins...
        #### move everything to place atom2 at the origin
        v3 = self.get_vector(atom3, check=False)
        self.translate_molecule(-v3)

        #### perform the actual rotation
        rot_matrix = compute_rotation_matrix(-self.get_vector(atom2, check=False), delta)

        for atom in atoms_to_move:
            self.geometry[atom] = np.dot(rot_matrix, self.get_vector(atom, check=False))

        error = self.get_dihedral(atom1, atom2, atom3, atom4, check=False) - dihedral

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
        for atom in range(1, self.num_atoms() + 1):
            self.geometry[atom] = self.geometry[atom] + vector

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

    def _recurse_through_formula(self, formula, masses, weights, cutoff=0.0000001, mass_precision=4, weight_precision=8):
        """
        Recurses through a formula and generates m/z isotopic pattern using tail recursion.

        To prevent blowup of memory, fragments with very low abundance are ignored. Masses and weights are also rounded after every step.
        To prevent error accumulation, internal precisions several orders of magnitude lower than the precision of interest should be employed.
        The default values should work nicely for low-res MS applications.

        Args:
            formula (np.array, dtype=np.int8): vector containing atoms left to incorporate
            masses (np.array): list of mass fragments at current iteration
            weights (np.array): relative weights at current iteration
            cutoff (float): cutoff for similarity (masses within ``cutoff`` will be combined)
            mass_precision (int): number of decimal places to store for mass
            weight_precision (int): number of decimal places to store for weight

        Returns:
            masses
            weights
        """
        if np.array_equal(formula, np.zeros(shape=92, dtype=np.int8)):
            return masses[np.argsort(masses)], weights[np.argsort(masses)]

        current_e = np.nonzero(formula)[0][0]
        e_masses, e_weights = get_isotopic_distribution(current_e)

        new_masses = np.zeros(shape=(len(masses)*len(e_masses)))
        new_weights = np.zeros(shape=(len(masses)*len(e_masses)))
        for i in range(len(masses)):
            for j in range(len(e_masses)):
                new_masses[i*len(e_masses)+j] = masses[i] + e_masses[j]
                new_weights[i*len(e_masses)+j] = weights[i] * e_weights[j]

        newer_masses, indices = np.unique(np.round(new_masses, decimals=mass_precision), return_inverse=True)
        newer_weights = np.zeros_like(newer_masses)
        for k in range(len(newer_weights)):
            newer_weights[k] = np.sum(new_weights[np.nonzero(indices == k)])
        newer_weights = np.round(newer_weights, decimals=weight_precision)

        formula[current_e] += -1
        above_cutoff = np.nonzero(newer_weights > cutoff)
        return self._recurse_through_formula(formula, newer_masses[above_cutoff], newer_weights[above_cutoff], cutoff, mass_precision, weight_precision)

    def calculate_mass_spectrum(self, **kwargs):
        """
        Generates list of m/z values based on formula string (e.g. "C10H12")

        Final weights rounded to one decimal point (because of low-res MS).
        """
        form_vec = np.zeros(shape=92, dtype=np.int8)
        for z in self.atomic_numbers:
            form_vec[z - 1] += 1

        masses, weights = self._recurse_through_formula(form_vec, [0], [1], **kwargs)

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

        if (not isinstance(coordinates, list)) or (len(coordinates) != 3):
            raise TypeError("coordinates must be list with three elements")

        if not isinstance(symbol, str):
            raise TypeError(f"symbol {symbol} must be a string!")

        number = get_number(symbol)
        self.atomic_numbers = np.append(self.atomic_numbers, [number]).view(OneIndexedArray)
        self.geometry = np.append(self.geometry, [coordinates], axis=0).view(OneIndexedArray)
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
            self.geometry = np.delete(self.geometry, number, axis=0).view(OneIndexedArray)
            self.atomic_numbers = np.delete(self.atomic_numbers, number).view(OneIndexedArray)
            return self
        except:
            raise ValueError("removing atom {number} failed!")

    def get_atomic_number(self, atom):
        """
        Get the atomic number for a given atom.

        Args:
            atom1 (int): number of the first atom

        Returns:
            the atomic number of that atom
        """
        self._check_atom_number(atom)
        return self.atomic_numbers[atom]

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

    def get_distance(self, atom1, atom2, check=True, _dist=compute_distance_between):
        """
        Wrapper to compute distance between two atoms.

        This function is relatively slow (rate-limiting for certain applications), so performance boosts have been implemented (e.g. preloading ``_dist``).

        Args:
            atom1 (int): number of the first atom
            atom2 (int): number of the second atom
            check (Bool): whether to validate input data (can be overridden to prevent slow double-checking)
            _dist (function): function usd to compute distance

        Returns:
            the distance, in Angstroms
        """
        if check:
            try:
                atom1 = int(atom1)
                atom2 = int(atom2)
            except:
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
            except:
                raise TypeError("atom numbers cannot be cast to int!")

            self._check_atom_number(atom1)
            self._check_atom_number(atom2)

        return np.sum(np.square(self.get_vector(atom1, atom2, check=False)))

    def get_angle(self, atom1, atom2, atom3, check=True, _angle=compute_angle_between):
        """
        Wrapper to compute angle between three atoms.

        This function is relatively slow (rate-limiting for certain applications), so performance boosts have been implemented (e.g. preloading ``_angle``).

        Args:
            atom1 (int): number of the first atom
            atom2 (int): number of the second atom
            atom3 (int): number of the third atom
            check (Bool): whether to validate input data (can be overridden to prevent slow double-checking)
            _angle (function): function usd to compute angle

        Returns:
            the angle, in degrees
        """
        if check:
            try:
                atom1 = int(atom1)
                atom2 = int(atom2)
                atom3 = int(atom3)
            except:
                raise TypeError("atom numbers cannot be cast to int!")

            self._check_atom_number(atom1)
            self._check_atom_number(atom2)
            self._check_atom_number(atom3)

        v1 = self.get_vector(atom1, check=False)
        v2 = self.get_vector(atom2, check=False)
        v3 = self.get_vector(atom3, check=False)

        return _angle(v1 - v2, v3 - v2)

    def get_dihedral(self, atom1, atom2, atom3, atom4, check=True, _dihedral=compute_dihedral_between):
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

        Returns:
            the dihedral angle, in degrees
        """
        if check:
            try:
                atom1 = int(atom1)
                atom2 = int(atom2)
                atom3 = int(atom3)
                atom4 = int(atom4)
            except:
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
        except:
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

    def optimize_dihedral(self, atom1, atom2, atom3, atom4):
        """
        Minimizes the value of ``self.rms_distance_between_atoms`` for the given dihedral, in one-degree increments.
        A cheap alternative to geometry optimization using *ab initio* methods or density functional theory.

        Args:
            atom1 (int): atom number of the first atom in the dihedral
            atom2 (int): atom number of the second atom in the dihedral
            atom3 (int): atom number of the third atom in the dihedral
            atom4 (int): atom number of the fourth atom in the dihedral

        Returns:
            the final value of the angle
        """
        self._check_atom_number(atom1)
        self._check_atom_number(atom2)
        self._check_atom_number(atom3)
        self._check_atom_number(atom4)

        best_angle = 0
        best_dist = 0

        for angle in range(0, 360):
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
        except:
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

