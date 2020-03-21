# import sys
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
        molecules (np.array): list of `Molecule` objects
        energies (np.array): list of energies
    """

    def __init__(self, name=None, **kwargs):
        """
        Create new instance, and optionally create a bunch of molecules too.

        Args:
            name (str): name of Ensemble
            **kwargs: to pass to ``self.batch_add()``
        """
        self.name = name
        self.molecules = np.array([])
        self.energies = np.array([])

        if all(arg in kwargs for arg in ["atomic_numbers", "geometries"]):
            self.batch_add(**kwargs)

    def __getitem__(self, key):
        return self.molecules[key]

    def __setitem__(self, key, item):
        self.molecules[key] = item

    def batch_add(self, atomic_numbers, geometries, bonds=None, charges=None, multiplicities=None):
        """
        Automatically generates ``Molecule`` objects and adds them using ``self.add_molecule()``.

        Args:
            atomic_numbers (list): list of lists of atomic numbers.
            geometry (list): list of 3-tuples of xyz coordinates
            bonds (list): list of edges (i.e. an n x 2 ``numpy`` array).
            charges (int): list of molecular charges - will default to 0 for all molecules.
            multiplicities (int): list of multiplicities - will default to 1 (singlet) for all molecules.
        """
        if (atomic_numbers is None) or (geometries is None):
            return

        if charges is None:
            charges = list(np.zeros(shape=len(atomic_numbers)))

        if multiplicities is None:
            multiplicities = list(np.ones(shape=len(atomic_numbers)))

        if bonds is None:
            bonds = [None] * len(atomic_numbers)

        assert all(
            len(x) == len(atomic_numbers) for x in [geometries, bonds, charges, multiplicities]
        ), "uneven list lengths -- cannot batch create molecules!"

        for numbers, geometry, bond_edges, charge, multiplicity in zip(atomic_numbers, geometries, bonds, charges, multiplicities):
            if len(numbers) != len(geometry):
                raise TypeError("atoms and geometries must be the same length!")

            mol = Molecule(numbers, geometry, bonds=bond_edges, charge=charge, multiplicity=multiplicity)
            self.add_molecule(mol)

    def add_molecule(self, molecule, energy=None):
        """
        Adds a molecule to the ensemble. ``copy.deepcopy`` is used so that an independent copy of the molecule is saved.

        Args:
            molecule (Molecule): the molecule to be added
        """
        if not isinstance(molecule, Molecule):
            raise TypeError("molecule is not a Molecule - so it can't be added!")

        self.molecules = np.append(self.molecules, [copy.deepcopy(molecule)])
        self.energies = np.append(self.energies, [energy])

    def _check_molecule_number(self, number):
        """
        Helper method which performs quick checks on the validity of a given molecule number.
        """
        try:
            number = int(number)
        except:
            raise TypeError(f"atom number {number} must be integer")

        if number >= len(self.molecules):
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
        use_energies = True
        for ensemble in ensembles:
            assert isinstance(ensemble, Ensemble), "can't join an object that isn't an Ensemble!"
            if len(ensemble.energies) != len(ensemble.molecules):
                use_energies = False

        for ensemble in ensembles:
            for idx, mol in np.ndenumerate(ensemble.molecules):
                if use_energies:
                    new_ensemble.add_molecule(mol, energy=ensemble.energies[idx])
                else:
                    new_ensemble.add_molecule(mol)

        return new_ensemble


class ConformationalEnsemble(Ensemble):
    """
    Class that represents a group of conformers. All members must have the same atom types in the same order.

    Allows you to align and remove redundant molecules, unlike ``Ensemble``.

    Attributes:
        name (str): name, for identification
        molecules (list): list of `Molecule` objects
    """

    def batch_add(self, atomic_numbers=None, geometries=None, **kwargs):
        """
        Automatically generates ``Molecule`` objects and adds them.

        Takes only a single molecule's value for ``charge``/``multiplicity``/``atomic_numbers``/``bonds``, not a list of lists like ``Ensemble.__init__()``.

        Args:
            atomic_numbers (list): list of atomic numbers.
            geometry (np.ndarray): list of 3-tuples of xyz coordinates
            bonds (list): list of edges (i.e. an n x 2 `numpy` array).
            charges (int): molecular charge - will default to 0 for all molecules.
            multiplicities (int): spin multiplicity - will default to 1 (singlet) for all molecules.
        """
        if (atomic_numbers is None) or (geometries is None):
            return

        for geometry, z in zip(geometries, atomic_numbers):
            if len(z) != len(geometry):
                raise TypeError("atoms and geometries must be the same length!")

        assert all(set(z) == set(atomic_numbers[0]) for z in atomic_numbers), "not all atomic numbers match; can't make ConformationalEnsemble"
        assert all(len(z) == len(atomic_numbers[0]) for z in atomic_numbers), "not all atomic numbers match; can't make ConformationalEnsemble"
        assert all(g.shape == geometries[0].shape for g in geometries), "not all geometries match; can't make ConformationalEnsemble!"

        for x in ["charges", "multiplicities", "bonds"]:
            if x in kwargs.keys():
                if isinstance(kwargs[x], list) and len(kwargs[x]) > 1:
                    assert all(y == kwargs[x][0] for y in kwargs[x]), f"not all {x} match; can't make ConformationalEnsemble!"

        super().batch_add(atomic_numbers, geometries, **kwargs)

    def add_molecule(self, molecule, energy=None):
        """
        Checks that the molecule contains the same atom types in the same order as existing molecules, and that the molecule has the same charge/multiplicity.
        """
        if len(self.molecules) > 0:
            if molecule.num_atoms() != self.molecules[0].num_atoms():
                raise ValueError("wrong number of atoms for this ensemble")

            if molecule.charge != self.molecules[0].charge:
                raise ValueError("wrong charge for this ensemble")

            if molecule.multiplicity != self.molecules[0].multiplicity:
                raise ValueError("wrong spin multiplicity for this ensemble")

            if not np.array_equal(molecule.atomic_numbers, self.molecules[0].atomic_numbers):
                raise ValueError("wrong atom types for this ensemble")

        super().add_molecule(molecule, energy=energy)

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
        use_energies = True
        for ensemble in ensembles:
            assert isinstance(ensemble, ConformationalEnsemble), "can't join an object that isn't an ConformationalEnsemble!"
            if len(ensemble.energies) != len(ensemble.molecules):
                use_energies = False

        for ensemble in ensembles:
            for idx, mol in np.ndenumerate(ensemble.molecules):
                if use_energies:
                    new_ensemble.add_molecule(mol, energy=ensemble.energies[idx])
                else:
                    new_ensemble.add_molecule(mol)

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
