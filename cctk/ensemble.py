import sys
import re
import numpy as np
import copy

from cctk import Molecule
from cctk.helper_functions import align_matrices, compute_RMSD


class Ensemble:
    """
    Class that represents a group of molecules. They do not all need to have the same atoms or bonds.

    Calling ``ensemble[x]`` is shorthand for calling ``ensemble.molecules[x]``.

    Attributes:
        name (str): name, for identification
        _items (dict):
            keys: ``Molecule`` objects
            values: dictionaries containing properties from each molecule, variable. should always be one layer deep.
    """

    def __init__(self, name=None):
        """
        Create new instance.

        Args:
            name (str): name of Ensemble
        """
        self.name = name
        self._items = {}

    def __getitem__(self, key):
        if isinstance(key, Molecule):
            return self._items[key]
        elif isinstance(key, int):
            return list(self._items)[key]
        else:
            raise KeyError(f"not a valid datatype for Ensemble key: {type(key)}")

    def __setitem__(self, key, item):
        if isinstance(key, Molecule):
            self._items[key] = item
        elif isinstance(key, int):
            list(self._items)[key] = item
        else:
            raise KeyError(f"not a valid datatype for Ensemble key: {type(key)}")

    def __len__(self):
        return len(self._items)

    def has_property(self, idx, prop):
        if prop in list(self._items[self[idx]].keys()):
            return True
        else:
            return False

    def items(self):
        return self._items.items()

    def add_molecule(self, molecule, properties={}):
        """
        Adds a molecule to the ensemble. ``copy.deepcopy`` is used so that an independent copy of the molecule is saved.

        Args:
            molecule (Molecule): the molecule to be added
        """
        if not isinstance(molecule, Molecule):
            raise TypeError("molecule is not a Molecule - so it can't be added!")

        mol = copy.deepcopy(molecule)
        self._items[mol] = properties

    def _check_molecule_number(self, number):
        """
        Helper method which performs quick checks on the validity of a given molecule number.
        """
        try:
            number = int(number)
        except:
            raise TypeError(f"atom number {number} must be integer")

        if number >= len(self._items):
            raise ValueError(f"atom number {number} too large!")

    @classmethod
    def join_ensembles(cls, ensembles, name=None):
        """
        Creates a new Ensemble object from existing ensembles.

        If every ensemble has energies defined, then the new ensemble will have energies defined too.

        Args:
            name (str): name of Ensemble created
            ensembles (list of Ensembles): Ensemble objects to join
        """
        new_ensemble = Ensemble(name=name)
        for ensemble in ensembles:
            assert isinstance(ensemble, Ensemble), "can't join an object that isn't an Ensemble!"

        for ensemble in ensembles:
            new_ensemble._items.update(ensemble.items)

        return new_ensemble

    def to_df(self):
        pass

class ConformationalEnsemble(Ensemble):
    """
    Class that represents a group of conformers. All members must have the same atom types in the same order.

    Allows you to align and remove redundant molecules, unlike ``Ensemble``.

    Attributes:
        name (str): name, for identification
        molecules (list): list of `Molecule` objects
    """

    def add_molecule(self, molecule, properties={}):
        """
        Checks that the molecule contains the same atom types in the same order as existing molecules, and that the molecule has the same charge/multiplicity.
        """
        if len(self._items) > 0:
            if molecule.num_atoms() != self[0].num_atoms():
                raise ValueError("wrong number of atoms for this ensemble")

            if molecule.charge != self[0].charge:
                raise ValueError("wrong charge for this ensemble")

            if molecule.multiplicity != self[0].multiplicity:
                raise ValueError("wrong spin multiplicity for this ensemble")

            if not np.array_equal(molecule.atomic_numbers, self[0].atomic_numbers):
                raise ValueError("wrong atom types for this ensemble")

        super().add_molecule(molecule, properties)

    @classmethod
    def join_ensembles(cls, ensembles, name=None):
        """
        Creates a new ConformationalEnsemble object from existing ensembles.

        If every ensemble has energies defined, then the new ensemble will have energies defined too.

        Args:
            name (str): name of ConformationalEnsemble created
            ensembles (list of ConformationalEnsembles): ConformationalEnsemble objects to join
        """
        new_ensemble = ConformationalEnsemble(name=name)
        for ensemble in ensembles:
            assert isinstance(ensemble, ConformationalEnsemble), "can't join an object that isn't an ConformationalEnsemble!"

        for ensemble in ensembles:
            for mol, prop in ensemble.items():
                    new_ensemble.add_molecule(mol, prop)

        return new_ensemble

    def align(self, align_to=0, atoms=None, return_rmsd=False):
        """
        Aligns every geometry to the specified geometry based on the atoms in `atom_numbers`. If `atom_numbers` is `None`, then a full alignment is performed.

        Args:
            align_to (int): which geometry to align to (0-indexed)
            atoms (list): which atoms to align in each molecule (1-indexed; must be at least 3)
                alternatively, specify ``None`` for all atoms or "heavy" for all heavy atoms
            return_rmsd (Bool): whether to return RMSD before and after rotation

        Returns:
            a new ``Ensemble()`` object with the objects aligned
            (optional) before rmsd and after rmsd
        """
        self._check_molecule_number(align_to)

        if atoms is None:
            atoms = np.arange(1, self.molecules[0].num_atoms() + 1)
        elif isinstance(atoms, str) and (atoms == "heavy"):
            atoms = self.molecules[0].get_heavy_atoms()
        else:
            try:
                atoms = np.array(atoms)
                if len(atoms) < 3:
                    raise ValueError("not enough atoms for alignment - need 3 in 3D space!")

            except:
                raise ValueError("atom numbers is not a recognized keyword and cannot be cast to numpy array... try again!")

        #### move everything to the center!
        for molecule in self.molecules:
            molecule.center()

        template = self.molecules[align_to].geometry[atoms]
        before_rmsd = 0
        after_rmsd = 0

        #### perform alignment using Kabsch algorithm
        new_ensemble = copy.deepcopy(self)
        for molecule in new_ensemble.molecules:
            before_rmsd += compute_RMSD(template, molecule.geometry[atoms])
            new_geometry = align_matrices(molecule.geometry[atoms], molecule.geometry, template)
            molecule.geometry = new_geometry
            after_rmsd += compute_RMSD(template, molecule.geometry[atoms])

            assert len(molecule.geometry) == len(molecule.atomic_numbers), "wrong number of geometry elements!"

        if return_rmsd:
            return new_ensemble, before_rmsd, after_rmsd
        else:
            return new_ensemble

    def eliminate_redundant(self, cutoff=0.5, heavy_only=True, atom_numbers=None):
        """
        Returns non-redundant conformations. When redundancies are found, only the first geometry is kept.
        This will change the numbering of all the ensembles!

        Args:
            cutoff (float): molecules with less than this value for RMSD will be considered redundant and eliminated.
            heavy_only (Bool): if ``True``, then only heavy atoms are considered for the RMSD calculation
            atom_numbers (list): 1-indexed list of atoms to consider for RMSD calculation - if present, overrides ``heavy_only``

        Returns:
            a new ``ConformationalEnsemble`` object where redundant conformers have been deleted and all molecules have been aligned
        """
        if atom_numbers:
            atom_numbers = [n - 1 for n in atom_numbers]
        else:
            if heavy_only:
                atom_numbers = self.molecules[0].get_heavy_atoms()
            else:
                atom_numbers = list(range(len(self.molecules[0])))

        for m in self.molecules:
            for n in atom_numbers:
                try:
                    #### atom_numbers is 0-indexed
                    m._check_atom_number(n + 1)
                except:
                    raise ValueError(f"molecule in ensemble does not have atom {n}!")

        #### align all molecules
        new_ensemble = self.align(atoms=atom_numbers)
        to_delete = [False] * len(new_ensemble.molecules)

        for i in range(len(new_ensemble.molecules)):
            if to_delete[i]:
                continue
            for j in range(i + 1, len(new_ensemble.molecules)):
                if to_delete[j]:
                    continue

                geometry1 = new_ensemble.molecules[i].geometry[atom_numbers]
                geometry2 = new_ensemble.molecules[j].geometry[atom_numbers]

                rmsd = compute_RMSD(geometry1, geometry2)
                if rmsd < cutoff:
                    to_delete[j] = True

        #### you have to delete in reverse order or you'll throw off the subsequent indices
        for i in sorted(range(len(new_ensemble.molecules)), reverse=True):
            if to_delete[i]:
                new_ensemble.molecules = np.delete(new_ensemble.molecules, i)
                new_ensemble.energies = np.delete(new_ensemble.energies, i)

        return new_ensemble

    def get_geometric_parameters(self, parameter, atom1, atom2, atom3=None, atom4=None):
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

    def get_lowest_energy(self, num=10):
        return self.molecules[np.argsort(self.energies)][0:10]

    def get_within_cutoff(self, cutoff=5):
        return self.molecules[self.energies <= (np.min(self.energies) + 5)]
