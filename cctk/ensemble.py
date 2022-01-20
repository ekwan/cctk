import numpy as np
from copy import deepcopy

import cctk
from cctk.helper_functions import align_matrices


class Ensemble:
    """
    Class representing a collection of molecules. They do not all need to have the same atoms or bonds.

    Ensembles are composed of molecules and properties. Molecules are ``Molecule`` objects, whereas properties are ``dict`` objects containing calculation-specific information.

    There are various shortcuts for handling ``Ensemble`` objects:

    - ``ensemble[molecule]`` or ``ensemble[0]`` will return new ``Ensemble`` objects with only the specified molecules.
        Lists or slices can also be used: so ``ensemble[0:10:2]`` or ``ensemble[[molecule1, molecule2, molecule3]]`` will also return new ``Ensemble`` objects.
    - Individual properties can be read through tuple indexing: ``ensemble[0,"energy"]`` will return the energy of the first molecule,
        while ``ensemble[:,"energy"]`` will return a list of all the energies.
    - To access ``Molecule`` objects, use ``ensemble.molecule``: ``ensemble.molecule[0]`` will return the first object, whereas ``ensemble.molecule[1:3]`` will return a list.
    - ``ensemble.items()`` will return a list of (molecule, property) pairs.
    - ``ensemble.molecule_list()`` and ``ensemble.properties_list()`` return lists of molecules and properties, respectively.

    Attributes:
        name (str): name, for identification
        _items (dict): keys: ``Molecule`` objects; values: dictionaries containing properties from each molecule, variable. should always be one layer deep.
        molecules (``MoleculeIndexer``): special object that accesses the keys
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
        return f"Ensemble (name={name}, {len(self._items)} molecules)"

    def __getitem__(self, key):
        if isinstance(key, (int, np.integer)):
            mol = self.molecule_list()[key]
            prop = self.properties_list()[key]
            new = type(self)(name=self.name) # will return either Ensemble or subclass thereof
            new.add_molecule(mol, properties=prop)
            return new
        elif isinstance(key, cctk.Molecule):
            idx = self.molecule_list().index(key)
            return self[idx]
        elif isinstance(key, (list, np.ndarray)):
            new_list = [self[k] for k in key]
            return self.join_ensembles(new_list, name=self.name)
        elif isinstance(key, slice):
            start, stop, step = key.indices(len(self))
            return self[list(range(start, stop, step))]
        elif isinstance(key, tuple):
            return self.get_property(key[0], key[1])
        elif key is None:
            return self
        else:
            raise KeyError(f"not a valid datatype for Ensemble key: {type(key)}")

    def __setitem__(self, key, item):
        assert isinstance(key, tuple), "need two indexes to set a value in an ensemble!"
        idx = key[0]
        name = key[1]

        if isinstance(idx, slice):
            start, stop, step = idx.indices(len(self))
            self[list(range(start, stop, step)), name] = item
        elif isinstance(idx, (list, np.ndarray)) and isinstance(item, (list, np.ndarray)):
            assert len(idx) == len(item), f"can't set {len(item)} items into {len(key)} variables (cf. pigeonhole principle)"
            for (k, i) in zip(idx, item):
                self[k, name] = i
        elif isinstance(idx, (list, np.ndarray)):
            for k in idx:
                self[k, name] = item
        elif isinstance(idx, (int, np.integer)):
            mol = self.molecule_list()[idx]
            self[mol, name] = item
        elif isinstance(idx, cctk.Molecule):
            if isinstance(name, (list, np.ndarray)):
                for n in name:
                    self[idx,n] = item
            #### we can't assign multiple items to a list of names since that would preclude assigning a list to a single variable
            else:
                self._items[idx][name] = item
        else:
            raise KeyError(f"not a valid datatype for Ensemble index: {type(idx)}")

    def __len__(self):
        return len(self._items)

    def __iter__(self):
        return iter(self.items())

    def keys(self):
        return self._items.keys()

    def values(self):
        return self._items.values()

    def molecule_list(self):
        """
        Returns a list of the constituent molecules.
        """
        return list(self.keys())

    def properties_list(self):
        """
        Returns a list of the constituent molecules.
        """
        return list(self.values())

    def has_property(self, idx, prop):
        """
        Returns ``True`` if property is defined for index ``idx`` and ``False`` otherwise.
        """
        combined = self.combined_properties()
        if prop in combined:
            return True
        else:
            return False

    def combined_properties(self):
        """
        Returns a dictionary containing the most up-to-date version of each property.
        """
        combined = dict()
        for p in self.properties_list():
            combined = {**combined, **p}
        return combined

    def get_property(self, idx, prop):
        """
        """
        ensemble = self[idx]
        result = []
        for m, p in ensemble.items():
            if isinstance(prop, list):
                row = []
                for x in prop:
                    if x in p:
                        row.append(p[x])
                    else:
                        row.append(None)
                result.append(row)
            else:
                if prop in p:
                    result.append(p[prop])
                else:
                    result.append(None)
        if len(ensemble) == 1:
            if result[0] is None:
                return None
            return result[0]
        else:
            found_something = False
            for x in result:
                if x is not None:
                    found_something = True
                    break
            if found_something:
                return result
            else:
                return None

    def get_properties_dict(self, idx):
        """
            Returns the dictionary of molecule properties for the specified molecule.

            Args:
                idx (int or cctk.Molecule): a molecule belonging to this ensemble, either
                                            0-indexed or given explicitly as a Molecule

            Returns:
                the property dict corresponding to this Molecule
        """
        assert isinstance(idx, (int, np.integer, cctk.Molecule)), "index must be int or Molecule"
        ensemble = self[idx]
        assert len(ensemble) == 1, "idx returned too many ensembles"
        return ensemble.properties_list()[0]

    def items(self):
        """
        Returns a list of (molecule, properties) tuple pairs.
        """
        return self._items.items()

    # object to allow convenient indexing of the molecules in the ensemble  
    #
    # allowed use cases
    #
    # retrieving molecules:
    # ensemble.molecules[0]: first molecule
    # ensemble.molecules[-1]: last molecule
    # ensemble.molecules[[0,1]]: first two molecules as a list
    # ensemble.molecules[0:4:2]: first and third molecules as a list
    #
    # setting molecule properties this way is not allowed
    class _MoleculeIndexer():
        def __init__(self, ensemble):
            self.ensemble = ensemble

        def __getitem__(self, key):
            items_list = list(self.ensemble._items.keys())
            n_items = len(items_list)
            if isinstance(key, (int, np.integer)):
                self._check_key(key, n_items)
                return items_list[key]
            if isinstance(key, np.ndarray):
                assert len(np.shape(key)) == 1, f"multidimensional keys not allowed, shape was {np.shape(key)}"
            if isinstance(key, (list, np.ndarray)):
                return_list = []
                for k in key:
                    assert isinstance(k, (int, np.integer)), f"key {k} in {str(key)} is not an integer, type is {str(type(k))}"
                    self._check_key(k, n_items)
                    return_list.append(items_list[k])
                return return_list
            elif isinstance(key, slice):
                start, stop, step = key.indices(n_items)
                return [ items_list[i] for i in range(start, stop, step) ]
            else:
                raise ValueError(f"cannot index with type {str(type(key))}")

        def __setitem__(self, key, item):
            raise LookupError("cannot set molecule properties this way; use ensemble.set_property_dict(molecule, property_dict) instead")

        def _check_key(self, key, n_items):
            assert -n_items <= key < n_items, f"key {key} is out of range...must be between {-n_items} and {n_items-1} inclusive"

        def __iter__(self):
            return iter(self.ensemble.molecule_list())

    def properties(self, num=None):
        """
        Returns a list of the constituent properties.
        """
        if num is None:
            return list(self.values())
        else:
            assert isinstance(num, int), "num must be integer"
            return list(self.values())[num]

    def sort_by(self, property_name, ascending=True):
        """
        Sorts the ensemble by the specified property.
        Throws an error if the property is missing for any entries.
        Consistent, sort-compatible property values are assumed and not checked.

        Args:
            property_name (str): the name of the property to sort on (must be a string or number)
            ascending (bool): whether the property should increase or decrease in value
        Returns:
            new Ensemble (current ensemble is not modified)
        """
        property_list = self[:,property_name]
        if property_list is None:
            raise ValueError(f"property '{property_name}' not found in ensemble")
        property_list = np.asarray(property_list)
        n_missing_entries = np.count_nonzero(property_list==None)
        if n_missing_entries > 0:
            error = "---sorting error---\n"
            error += str(property_list)
            raise ValueError(f"{error}\nproperty '{property_name}' has {n_missing_entries} missing entries and cannot be sorted")
        new_indices = np.argsort(property_list)
        if not ascending:
            new_indices = np.flip(new_indices)
        return self[[new_indices]]

    def add_molecule(self, molecule, properties=None, copy=False):
        """
        Adds a molecule to the ensemble.

        Args:
            molecule (Molecule): the molecule to be added
            properties (dict): property name (str) to property value
            copy (bool): whether to store an independent copy of the molecule
        """
        if not isinstance(molecule, cctk.Molecule):
            raise TypeError("molecule is not a Molecule - so it can't be added!")

        if copy:
            molecule = deepcopy(molecule)

        if properties is None:
            #### empty dicts all point to the same memory address by default, so need to prevent that behavior by initializing non-empty dict
            properties = {"placeholder": 1}
            del properties["placeholder"]

        assert isinstance(properties, dict), f"properties must be a dict and not type {type(properties)}"

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
            new_ensemble._items.update(ensemble.items())

        return new_ensemble

    def lowest_molecules(self, property_name, num=1):
        """
        Retrieves the molecules with the lowest values of the specified property.

        Args:
            property_name (str): the name of the property to sort on
            num (int): how many molecules to return
        Returns:
            lowest ``Molecule`` (if num==1)
            ``list`` of ``Molecule`` (otherwise)
        """
        assert isinstance(num, (int, np.integer)), f"num must be an integer, got {type(num)}"
        assert num > 0, f"num must be > 0, got {num}"
        sorted_ensemble = self.sort_by(property_name)
        if num > 1:
            return sorted_ensemble.molecules[0:num]
        return sorted_ensemble.molecules[0]

class ConformationalEnsemble(Ensemble):
    """
    Class that representing a group of conformers. All members must have the same atom types in the same order.
    """

    def __str__(self):
        n_atoms = 0
        if len(self._items) > 0:
            first_molecule = self.molecule_list()[0]
            n_atoms = first_molecule.num_atoms()
        if self.name is not None:
            return f"ConformationalEnsemble (name={self.name}, {len(self._items)} molecules, {n_atoms} atoms)"
        else:
            return f"ConformationalEnsemble ({len(self._items)} molecules, {n_atoms} atoms)"

    def add_molecule(self, molecule, properties=None, copy=False, checks=True):
        """
        Checks that the molecule contains the same atom types in the same order as existing molecules, and that the molecule has the same charge/multiplicity.
        """
        if len(self._items) > 0:
            initial_mol = self.molecule_list()[0]
            if molecule.num_atoms() != initial_mol.num_atoms():
                raise ValueError("wrong number of atoms for this ensemble")

            if molecule.charge != initial_mol.charge:
                raise ValueError("wrong charge for this ensemble")

            if molecule.multiplicity != initial_mol.multiplicity:
                raise ValueError("wrong spin multiplicity for this ensemble")

            if checks and not np.array_equal(molecule.atomic_numbers, initial_mol.atomic_numbers):
                raise ValueError("wrong atom types for this ensemble")

            #### only save one copy to save space
            molecule.bonds = initial_mol.bonds
            molecule.atomic_numbers = initial_mol.atomic_numbers

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
        n_atoms = self.molecules[0].num_atoms()

        if isinstance(comparison_atoms, str):
            if comparison_atoms == "all":
                comparison_atoms = np.arange(1, n_atoms + 1)
            elif comparison_atoms == "heavy":
                comparison_atoms = self.molecules[0].get_heavy_atoms()
        assert isinstance(comparison_atoms, (list, np.ndarray, cctk.OneIndexedArray)), f"unexpected type for comparison_atoms: {str(type(comparison_atoms))}"
        for a in comparison_atoms:
            assert 1 <= a <= n_atoms, f"atom number out of range: got {a}, but must be between 1 and {n_atoms}"

        assert len(comparison_atoms) >= 3, f"need at least 3 atoms for alignment, but only got {len(comparison_atoms)}"

        # duplicate the ensemble
        new_ensemble = deepcopy(self)

        # translate all molecules to the origin
        # with respect to the comparison atoms
        for molecule, _ in new_ensemble:
            full_geometry = molecule.geometry
            partial_geometry = full_geometry[comparison_atoms]
            translation_vector = -partial_geometry.mean(axis=0)
            molecule.translate_molecule(translation_vector)

        full_template_geometry = new_ensemble.molecules[to_geometry].geometry
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
                partial_geometry = new_geometry[comparison_atoms]
                after_RMSD = cctk.helper_functions.compute_RMSD(partial_template_geometry, partial_geometry)
                after_RMSDs.append(after_RMSD)
            assert len(molecule.geometry) == n_atoms, f"wrong number of geometry elements! expected {n_atoms}, got {len(molecule.geometry)}"

        if compute_RMSD:
            return new_ensemble, before_RMSDs, after_RMSDs
        return new_ensemble

    def eliminate_redundant(self, RMSD_cutoff=0.5, comparison_atoms="heavy", return_RMSD=False):
        """
        Aligns every geometry in this ensemble and then creates a new ensemble that contains only the non-redundant conformers.
        If energies are available, the lowest energy conformer will be kept for every redundancy.
        The current ensemble will not be modified.  The resulting ensemble will be sorted by energy (if available).

        Args:
            RMSD_cutoff (float): remove conformers that are more similar than this threshold
            to_geometry (int): the reference geometry to align to (0-indexed)
            comparison_atoms (str or list): which atoms to use when computing alignments
                                            "heavy" for all non-hydrogen atoms,
                                            "all" for all atoms, or
                                            a list of 1-indexed atom numbers
            return_RMSD (bool): whether or not to return list of RMSD values

        Returns:
            new ``ConformationalEnsemble``, RMSDs to the reference geometry
        """
        # check inputs
        n_atoms = self.molecules[0].num_atoms()
        if isinstance(comparison_atoms, str):
            if comparison_atoms == "all":
                comparison_atoms = np.arange(1, n_atoms + 1)
            elif comparison_atoms == "heavy":
                comparison_atoms = self.molecules[0].get_heavy_atoms()

        assert isinstance(comparison_atoms, (list, np.ndarray, cctk.OneIndexedArray)), f"unexpected type for comparison_atoms: {str(type(comparison_atoms))}"
        for a in comparison_atoms:
            assert 1 <= a <= n_atoms, f"atom number out of range: got {a}, but must be between 1 and {n_atoms}"
        assert len(comparison_atoms) >= 3, f"need at least 3 atoms for alignment, but only got {len(comparison_atoms)}"

        assert isinstance(RMSD_cutoff, (float, int)), f"RMSD cutoff must be a float but got {str(type(RMSD_cutoff))}"
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

        # boolean indexing noticeably faster
        idxs = np.array(comparison_atoms)
        mask = np.zeros(old_ensemble.molecules[0].geometry.shape[0], dtype=bool)
        mask[idxs - 1] = True

        partial_geoms = [m.geometry[mask] for m in old_ensemble.molecules]
        new_partial_geoms = []

        rmsds = list()

        # add molecules one by one
        new_ensemble = ConformationalEnsemble()
        for i in sorted_indices:
            ok_to_add = True

            candidate_rmsd = 0
            for existing_molecule in new_partial_geoms:
                candidate_rmsd = cctk.helper_functions.compute_RMSD(partial_geoms[i], existing_molecule, checks=False)
                if candidate_rmsd < RMSD_cutoff:
                    ok_to_add = False
                    break

            if ok_to_add:
                candidate_molecule = old_ensemble.molecules[i]
                candidate_molecule_properties = old_ensemble.get_properties_dict(candidate_molecule)

                new_ensemble.add_molecule(candidate_molecule, candidate_molecule_properties)
                new_partial_geoms.append(candidate_molecule.geometry[mask])
                rmsds.append(candidate_rmsd)

        if return_RMSD:
            return new_ensemble, rmsds
        else:
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
        output = [None] * len(self)
        for index, molecule in enumerate(self.molecule_list()):
            if parameter == "distance":
                output[index] = molecule.get_distance(atom1, atom2)
            elif parameter == "angle":
                if atom3 is None:
                    raise ValueError("need atom3 to calculate angle!")
                output[index] = molecule.get_angle(atom1, atom2, atom3)
            elif parameter == "dihedral":
                if (atom3 is None) or (atom4 is None):
                    raise ValueError("need atom3 and atom4 to calculate dihedral!")
                output[index] = molecule.get_dihedral(atom1, atom2, atom3, atom4)
            else:
                raise ValueError(f"Invalid parameter {parameter}!")

        return output

    def assign_connectivity(self, index=0):
        """
        Assigns connectivity for all molecules based on molecule of index ``index``. Much faster than assigning connectivity for each individually -- but assumes all bonding is the same.
        """
        assert isinstance(index, int), "Need integer index"
        bonds = self.molecules[index].assign_connectivity().bonds

        for mol in self.molecules:
            mol.bonds = bonds

        return self

    def boltzmann_average(self, which, energies=None, temp=298, energy_unit="hartree", return_weights=False):
        """
        Computes the Boltzmann-weighted average of a property over the whole ensemble.

        Args:
            which (str): which property to compute
            energy (np.ndarray): list of energies to use for weighting.
                Will default to ``self[:,"energy"]``, although other strings can be passed as well as shorthand for ``self[:,energy]``.
            temp (float): temperature for Boltzmann-weighting, in K
            energy_unit (str): either ``kcal_mol`` or ``hartree``
            return_weights (bool): whether to return a list of weights too

        Returns:
            weighted property, of the same shape as the individual property
        """
        if energies is None:
            energies = self[:,"energy"]
        elif isinstance(energies, str):
            energies = self[:,energies]
        elif isinstance(energies, (list, np.ndarray, cctk.OneIndexedArray)):
            pass
        else:
            raise ValueError(f"invalid energy value {energies} (type {type(energies)})")

        for i, (m, pd) in enumerate(self.items()):
            assert which in pd, f"molecule #{i} doesn't have property {which} defined!"

        values = np.array(self[:,which], dtype=np.float64)
        energies = np.array(energies, dtype=np.float64)

        assert len(energies) == len(self)
        assert len(values) == len(self)
        assert all([e is not None for e in energies]), "energy not defined for all molecules"
        assert all([v is not None for v in values]), f"property {which} not defined for all molecules"

        # perhaps at some point we will need a real unit system like simtk/OpenMM, but not today!
        if energy_unit == "kcal_mol":
            energies = energies / 627.509
        energies = energies - np.min(energies)

        R = 3.1668105e-6 # eH/K

        weights = np.exp(-1*energies/(R*temp))
        weights = weights / np.sum(weights)

        try:
            weighted_value = np.average(values, weights=weights)
        except Exception as e:
            raise ValueError(f"error computing Boltzmann average: {e}")

        if return_weights:
            return weighted_value, weights
        else:
            return weighted_value
