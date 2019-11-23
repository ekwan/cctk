import sys
import re
import numpy as np
import copy

from cctk import Molecule
from cctk.helper_functions import align_matrices, compute_RMSD

class Ensemble():
    """
    Class that represents a group of molecules.

    Attributes:
        name (str): name, for identification
        molecules (list): list of `Molecule` objects
    """

    def __init__(self, name=None, atoms=None, geometries=None, bonds=None, charge=0, multiplicity=1):
        self.name = name
        self.molecules = []

        if atoms:
            self.batch_add(atoms, geometries, bonds, charge, multiplicity)

    def batch_add(self, atoms, geometries, bonds, charge=0, multiplicity=1):
        """
        Automatically generates ``Molecule`` objects and adds them.

        Args:
            atoms (list): list of atomic symbols  (same for each ensemble member)
            geometry (list): Numpy array of 3-tuples of xyz coordinates
            bonds (list): list of edges (i.e. an n x 2 `numpy` array). Same for each ensemble member.
            charge (int): the charge of the molecule
            multiplicity (int): the spin state of the molecule (1 corresponds to singlet, 2 to doublet, 3 to triplet, etc. -- so a multiplicity of 1 is equivalent to S=0)
        """
        for geometry in geometries:
            print(atoms)
            print(geometry)
            if len(atoms) != len(geometry):
                raise TypeError("atoms and geometries must be the same length!")

            mol = Molecule(atoms, geometry, bonds=bonds, charge=charge, multiplicity=multiplicity)
            self.add_molecule(mol)

    def add_molecule(self, molecule):
        """
        Adds a molecule to the ensemble. `copy.deepcopy` is used so that an independent copy of the molecule is saved.

        Checks that the molecule contains the same atom types in the same order as existing molecules.

        Args:
            molecule (Molecule): the molecule to be added
        """
        if not isinstance(molecule, Molecule):
            raise TypeError("molecule is not a Molecule - so it can't be added!")

        if len(self.molecules) > 0:
            if len(molecule.atoms) != len(self.molecules[0].atoms):
                raise ValueError("wrong number of atoms for this ensemble")

            if not np.array_equal(molecule.atoms, self.molecules[0].atoms):
                raise ValueError("wrong atom types for this ensemble")

        self.molecules.append(copy.deepcopy(molecule))

    def align (self, align_to=1, atoms=None):
        """
        Aligns every geometry to the specified geometry based on the atoms in `atom_numbers`. If `atom_numbers` is `None`, then a full alignment is performed.

        Args:
            align_to (int): which geometry to align to (1-indexed)
            atoms (list): which atoms to align in each molecule (1-indexed; must be at least 3)
        """
        self._check_molecule_number(align_to)

        if atoms and (len(atoms) < 3):
            raise ValueError("not enough atoms for alignment - need 3 in 3D space!")

        try:
            atoms = np.array(atoms)
            atoms += -1
        except:
            raise ValueError("atom_numbers cannot be cast to numpy array... disappointing!")

        #### atom numbers is 0-indexed now
        #### move everything to the center!
        for molecule in self.molecules:
            centroid = molecule.geometry[atoms].mean(axis=0)
            molecule.translate_molecule(-centroid)

        template = self.molecules[align_to-1].geometry[atoms]

        #### perform alignment
        for molecule in self.molecules:
            new_geometry = align_matrices(molecule.geometry[atoms], molecule.geometry, template)
            self.geometry = new_geometry

            assert len(molecule.geometry) == len(molecule.atoms), "wrong number of geometry elements!"

    def eliminate_redundant(self, cutoff=0.5, heavy_only=True, atom_numbers=None):
        """
        Returns non-redundant conformations. When redundancies are found, only the first geometry is kept.
        This will change the numbering of all the ensembles!

        Args:
            cutoff (float): molecules with less than this value for RMSD will be considered redundant and eliminated.
            heavy_only (Bool): if `True`, then only heavy atoms are considered for the RMSD calculation
            atom_numbers (list): 1-indexed list of atoms to consider for RMSD calculation - if present, overrides `heavy_only`
        """
        if atom_numbers:
            atom_numbers = [n-1 for n in atom_numbers]
        else:
            if heavy_only:
                atom_numbers = self.molecules[0].get_heavy_atoms()
            else:
                atom_numbers = list(range(len(self.molecules[0])))

        for m in self.molecules:
            for n in atom_numbers:
                try:
                    #### atom_numbers is 0-indexed
                    m._check_atom_number(n+1)
                except:
                    raise ValueError(f"molecule in ensemble does not have atom {n}!")

        #### align all molecules 
        self.align(atoms=atom_numbers)

        to_delete = [False] * len(self.molecules)

        for i in range(len(self.molecules)):
            if to_delete[i]:
                continue
            for j in range(i+1, len(self.molecules)):
                if to_delete[j]:
                    continue

                geometry1 = self.molecules[i].geometry[atom_numbers]
                geometry2 = self.molecules[j].geometry[atom_numbers]

                rmsd = compute_RMSD(geometry1, geometry2)
                if rmsd < cutoff:
                    to_delete[j] = True

        #### you have to delete in reverse order or you'll throw off the subsequent indices 
        for i in sorted(range(len(self.molecules)), reverse=True):
            if to_delete[i]:
                del self.molecules[i]

    def _check_molecule_number(self, number):
        """
        Helper method which performs quick checks on the validity of a given molecule number.
        """
        try:
            number = int(number)
        except:
            raise TypeError(f"atom number {number} must be integer")

        if number > len(self.molecules):
            raise ValueError(f"atom number {number} too large!")

        if number <= 0:
            raise ValueError(f"atom number {number} invalid: must be a positive integer!")

    def print_geometric_parameters(self, parameter, atom1, atom2, atom3=None, atom4=None):
        """
        Computes and outputs geometric parameters (bond distances, angles, or dihedral angles) for every member of ``self.molecules.``

        Args:
            parameter (str): one of ``angle``, ``distance``, or ``dihedral``
            atom1 (int): number of the atom in question
            atom2 (int): same, but for the second atom
            atom3 (int): same, but for the third atom (only required for parameter ``angle`` or ``dihedral``)
            atom4 (int): same, but for the fourth atom (only required for parameter ``dihedral``)

        Returns:
            a list of the specified parameter's values for each geometry
        """
        output = [None] * len(self.molecules)
        for index, molecule in enumerate(self.molecules):
            if parameter == "distance":
                output[index] = molecule.get_distance(atom1, atom2)
            elif parameter == "angle":
                if atom3 == None:
                    raise ValueError("need atom3 to calculate angle!")
                output[index] = molecule.get_angle(atom1, atom2, atom3)
            elif parameter == "dihedral":
                if (atom3 == None) or (atom4 == None):
                    raise ValueError("need atom3 and atom4 to calculate dihedral!")
                output[index] = molecule.get_dihedral(atom1, atom2, atom3, atom4)
            else:
                ValueError("Invalid parameter {}!".format(parameter))

        return output
