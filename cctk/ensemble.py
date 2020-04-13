import sys
import re
import numpy as np
import copy

import cctk
from cctk import Molecule
from cctk.helper_functions import align_matrices


class Ensemble:
    """
    Class that represents a group of molecules. They do not all need to have the same atoms or bonds.

    Ensembles are composed of molecules and properties. Molecules are ``Molecule`` objects, whereas properties are ``dict`` objects containing calculation-specific information.

    There are various shortcuts for handling ``Ensemble`` objects:
    - ``Ensemble`` can be treated like a list of ``Molecule`` objects, so ``ensemble[0]`` will return the first molecule, ``len(ensemble)`` the number of molecules, and so forth.
    - ``Ensemble`` can also return the corresponding properties, when passed a molecule: so ``ensemble[molecule]`` will return the corresponding properties dictionary.
    - ``Ensemble`` can also take tuples, allowing for dataframe-like behavior: so ``ensemble[1:3, "energy"]`` will return the energies of molecules 2-4.

    Attributes:
        name (str): name, for identification
        _items (dict):
            keys: ``Molecule`` objects
            values: dictionaries containing properties from each molecule, variable. should always be one layer deep.
            molecules: special object that accesses the keys
    """

    def __init__(self, name=None):
        """
        Create new instance.

        Args:
            name (str): name of Ensemble
        """
        self.name = name
        self._items = {}
        self.molecules = self._MoleculeIndexer(self)

    def __str__(self):
        name = "None" if self.name is None else self.name
        return f"Ensemble (name={name}, {len(_items)} molecules)"

    def __getitem__(self, key):
        if isinstance(key, Molecule):
            return self._items[key]
        elif isinstance(key, (int,np.integer)):
            return list(self._items)[key]
        elif isinstance(key, list):
            return [self[k] for k in key]
        elif isinstance(key, slice):
            start, stop, step = key.indices(len(self))
            return [self[i] for i in range(start, stop, step)]
        elif isinstance(key, tuple):
            (key1, key2) = key
            mols = self[key1]
            if isinstance(key1, list):
                if all(isinstance(k, Molecule) for k in key1):
                    mols = key1
            elif isinstance(key1, Molecule):
                mols = key1
            if key2 is None:
                return mols
            else:
                try:
                    if isinstance(mols, list):
                        return [self[mol][key2] for mol in mols]
                    else:
                        return self[mols][key2]
                except KeyError as e:
                    raise KeyError(f"property {key2} not defined for all molecules!")
        else:
            raise KeyError(f"not a valid datatype for Ensemble key: {type(key)}")

    def __setitem__(self, key, item):
        if isinstance(key, Molecule):
            assert isinstance(item, dict), "properties must be dict"
            self._items[key] = item
        elif isinstance(key, int):
            assert isinstance(item, Molecule), "can't add something that isn't a molecule to keys of ensemble!"
            list(self._items)[key] = item
        elif isinstance(key, list):
            assert len(item) == len(key), "wrong number of items to assign to list index!"
            for k, i in zip(key, item):
                self[k] = i
        elif isinstance(key, slice):
            start, stop, step = key.indices(len(self))
            assert len(item) == len(range(start, stop, step)), "wrong number of items to assign to slice index!"
            for i in range(start, stop, step):
                self[i] = item[i]
        elif isinstance(key, tuple):
            (key1, key2) = key
            mols = self[key1]
            if isinstance(key1, list):
                if all(isinstance(k, Molecule) for k in key1):
                    mols = key1
            elif isinstance(key1, Molecule):
                mols = key1
            if isinstance(mols, list):
                assert len(item) == len(mols), "wrong number of items to assign to tuple index"
                for m, i in zip(mols, item):
                    self[m][key2] = i
            else:
                self[mols][key2] = item
        else:
            raise KeyError(f"not a valid datatype for Ensemble key: {type(key)}")

    def __len__(self):
        return len(self._items)

    def __iter__(self):
        return iter(self.items())

    def keys(self):
        return self._items.keys()

    def values(self):
        return self._items.values()

    def has_property(self, idx, prop):
        """
        Returns ``True`` if property is defined for index ``idx`` and ``False`` otherwise.
        """
        if prop in list(self._items[self[idx]].keys()):
            return True
        else:
            return False

    def items(self):
        """
        Returns a list of (molecule, properties) tuple pairs.
        """
        return self._items.items()

    #def molecules(self):
    #    """
    #    Returns a list of the constituent molecules.
    #    """
    #    return list(self.keys())

    class _MoleculeIndexer():
        def __init__(self, ensemble):
            self.ensemble = ensemble

        def __getitem__(self, key):
            items_list = list(self.ensemble._items.keys())
            n_items = len(items_list)
            if isinstance(key, int):
                return items_list[key]
            elif isinstance(key, list):
                return_list = []
                for k in key:
                    assert isinstance(k, (int, np.integer)), f"key {k} is not an integer, type is {str(type(k))}"
                    assert 0 <= k < n_items, f"key {k} is out of range...must be between 0 and {n_items-1}"
                    return_list.append(items_list[k])
                return return_list
            elif isinstance(key, slice):
                start, stop, step = key.indices(n_items)
                return [ items_list[i] for i in range(start, stop, step) ]
            else:
                raise ValueError(f"cannot index with type {str(type(key))}")

        def __setitem__(self, key, item={}):
            assert isinstance(item, dict), f"expected item to be a dict but got {str(type(item))}"
            items_list = list(self.ensemble._items.keys())
            n_items = len(items_list)
            if isinstance(key, int):
                m = items_list[key]
                self.ensemble._items[m] = item
            elif isinstance(key, Molecule):
                assert key in self.ensemble, "molecule not in ensemble"
                self.ensemble._items[key] = item
            else:
                raise ValueError(f"not a valid type for setting an item: {type(key)}")

    def properties(self):
        """
        Returns a list of the constituent properties.
        """
        return list(self.values())

    def add_molecule(self, molecule, properties={}, copy=False):
        """
        Adds a molecule to the ensemble.

        Args:
            molecule (Molecule): the molecule to be added
            properties (dict): property name (str) to property value
            copy (bool): whether to store an independent copy of the molecule
        """
        if not isinstance(molecule, Molecule):
            raise TypeError("molecule is not a Molecule - so it can't be added!")
        assert isinstance(properties, dict), "properties must be a dict"

        if copy:
            molecule = copy.deepcopy(molecule)
        self._items[molecule] = properties

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


class ConformationalEnsemble(Ensemble):
    """
    Class that represents a group of conformers. All members must have the same atom types in the same order.

    """

    def __str__(self):
        n_atoms = 0
        if len(self._items) > 0:
            first_molecule = self[0]
            n_atoms = first_molecule.num_atoms()
        if self.name is not None:
            return f"ConformationalEnsemble (name={self.name}, {len(self._items)} molecules, {n_atoms} atoms)"
        else:
            return f"ConformationalEnsemble ({len(self._items)} molecules, {n_atoms} atoms)"

    def add_molecule(self, molecule, properties={}, copy=False):
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

            #### only save one copy to save space
            molecule.bonds = self[0].bonds
            molecule.atomic_numbers = self[0].atomic_numbers

        super().add_molecule(molecule, properties, copy)

    @classmethod
    def join_ensembles(cls, ensembles, name=None, copy=False):
        """
        Creates a new ConformationalEnsemble object from existing ensembles.
        Both molecules and properties are copied.

        Args:
            name (str): name of ConformationalEnsemble created
            ensembles (list of ConformationalEnsembles): ConformationalEnsemble objects to join
            copy (bool): whether to make copies of the component molecules
        """
        new_ensemble = ConformationalEnsemble(name=name)
        for ensemble in ensembles:
            assert isinstance(ensemble, ConformationalEnsemble), "can't join an object that isn't a ConformationalEnsemble!"

        for ensemble in ensembles:
            for mol, prop in ensemble.items():
                    new_ensemble.add_molecule(mol, prop, copy)

        return new_ensemble

    def align(self, to_geometry=0, comparison_atoms="heavy", compute_RMSD=False):
        """
        Aligns every geometry in this ensemble to the specified geometry,
        optionally computing the root-mean-square distance between each
        geometry and the reference geometry.

	Alignments are based on `atom_numbers`.
        The current ensemble will not be altered.  RMSDs will be calculated over the
        comparison atoms only.

        Args:
            to_geometry (int): the reference geometry to align to (0-indexed)
            comparison_atoms (str or list): which atoms to use when computing alignments
                                            "heavy" for all non-hydrogen atoms,
                                            "all" for all atoms, or
                                            a list of 1-indexed atom numbers
            compute_RMSD (Bool): whether to return RMSD before and after rotation

        Returns:
            new aligned ``ConformationalEnsemble`` or
            new aligned ``ConformationalEnsemble``, before_RMSD array, after_RMSD array
        """
        # check inputs
        self._check_molecule_number(to_geometry)
        n_atoms = self[0].num_atoms()

        if isinstance(comparison_atoms, str):
            if comparison_atoms == "all":
                comparison_atoms = np.arange(1, n_atoms + 1)
            elif comparison_atoms == "heavy":
                comparison_atoms = self[0].get_heavy_atoms()
        assert isinstance(comparison_atoms, (list, np.ndarray, cctk.OneIndexedArray)), f"unexpected type for comparison_atoms: {str(type(comparison_atoms))}"
        for a in comparison_atoms:
            assert 1 <= a <= n_atoms, f"atom number out of range: got {a}, but must be between 1 and {n_atoms}"

        assert len(comparison_atoms) >= 3, f"need at least 3 atoms for alignment, but only got {len(comparison_atoms)}"

        # duplicate the ensemble
        new_ensemble = copy.deepcopy(self)

        # translate all molecules to the origin
        # with respect to the comparison atoms
        for molecule, _ in new_ensemble:
            full_geometry = molecule.geometry
            partial_geometry = full_geometry[comparison_atoms]
            translation_vector = -partial_geometry.mean(axis=0)
            molecule.translate_molecule(translation_vector)

        full_template_geometry = new_ensemble[to_geometry].geometry
        partial_template_geometry = full_template_geometry[comparison_atoms]
        before_RMSDs = []
        after_RMSDs = []

        # perform alignment using Kabsch algorithm
        for i, (molecule, _) in enumerate(new_ensemble):
            full_geometry = molecule.geometry
            partial_geometry = full_geometry[comparison_atoms]
            if compute_RMSD:
                before_RMSD = cctk.helper_functions.compute_RMSD(partial_template_geometry, partial_geometry)
                before_RMSDs.append(before_RMSD)
            new_geometry = align_matrices(partial_geometry, full_geometry, partial_template_geometry)
            molecule.geometry = new_geometry
            if compute_RMSD:
                after_RMSD = cctk.helper_functions.compute_RMSD(new_ensemble[to_geometry], new_ensemble[i], comparison_atoms)
                after_RMSDs.append(after_RMSD)
            assert len(molecule.geometry) == n_atoms, f"wrong number of geometry elements! expected {n_atoms}, got {len(molecule.geometry)}"

        if compute_RMSD:
            return new_ensemble, before_RMSDs, after_RMSDs
        return new_ensemble

    def eliminate_redundant(self, RMSD_cutoff=0.5, comparison_atoms="heavy"):
        """
        Aligns every geometry in this ensemble to the specified geometry,
        and then creates a new ensemble that contains only the non-redundant conformers.
        If energies are available, the lowest energy conformer will be kept for every redundancy.
        The current ensemble will not be modified.

        Args:
            RMSD_cutoff (float): remove conformers that are more similar than this threshold
            to_geometry (int): the reference geometry to align to (0-indexed)
            comparison_atoms (str or list): which atoms to use when computing alignments
                                            "heavy" for all non-hydrogen atoms,
                                            "all" for all atoms, or
                                            a list of 1-indexed atom numbers

        Returns:
            new ```ConformationalEnsemble```, RMSDs to the reference geometry
        """
        # check inputs
        n_atoms = self[0].num_atoms()
        if isinstance(comparison_atoms, str):
            if comparison_atoms == "all":
                comparison_atoms = np.arange(1, n_atoms + 1)
            elif comparison_atoms == "heavy":
                comparison_atoms = self[0].get_heavy_atoms()
        assert isinstance(comparison_atoms, (list, np.ndarray, cctk.OneIndexedArray)), f"unexpected type for comparison_atoms: {str(type(comparison_atoms))}"
        for a in comparison_atoms:
            assert 1 <= a <= n_atoms, f"atom number out of range: got {a}, but must be between 1 and {n_atoms}"

        assert len(comparison_atoms) >= 3, f"need at least 3 atoms for alignment, but only got {len(comparison_atoms)}"

        assert isinstance(RMSD_cutoff, float), f"RMSD cutoff must be a float but got {str(type(RMSD_cutoff))}"
        assert RMSD_cutoff > 0.0001, "must use a big enough RMSD cutoff"

        # align all molecules
        old_ensemble = self.align(to_geometry=0, comparison_atoms=comparison_atoms, compute_RMSD=False)

        # sort molecules by energy if available
        energies_available = True
        for molecule,properties in old_ensemble.items():
            if "energy" not in properties:
                energies_available = False
                break

        n_molecules = len(old_ensemble)
        sorted_indices = list(range(n_molecules))
        if energies_available:
            energies = old_ensemble[:,"energy"]
            sorted_indices = list(np.argsort(energies))

        # add molecules one by one
        new_ensemble = ConformationalEnsemble()
        for i in sorted_indices:
            candidate_molecule = old_ensemble[i]
            candidate_molecule_properties = old_ensemble[candidate_molecule]
            ok_to_add = True
            for existing_molecule in new_ensemble.molecules():
                candidate_rmsd = cctk.helper_functions.compute_RMSD(candidate_molecule, existing_molecule, comparison_atoms, checks=False)
                if candidate_rmsd < RMSD_cutoff:
                    ok_to_add = False
                    break
            if ok_to_add:
                new_ensemble.add_molecule(candidate_molecule, candidate_molecule_properties)
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
