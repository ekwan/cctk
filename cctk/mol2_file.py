import re
import numpy as np
import networkx as nx

from cctk import File, Ensemble, ConformationalEnsemble, Molecule
from cctk.helper_functions import get_symbol, get_number


class MOL2File(File):
    """
    Class representing SYBYL ``.mol2`` files.

    Attributes:
        name (str): name of file
        ensemble (Ensemble): ``Ensemble`` or ``ConformationalEnsemble`` object
    """

    def __init__(self, name=None):
        if isinstance(name, str):
            self.name = name

    @classmethod
    def read_file(cls, filename, name=None, **kwargs):
        """
        Reads ``.mol2`` file and generates a ``MOL2File`` instance.

        Args:
            filename (str): path to file
            name (str): name of the file

        Returns:
            MOL2File object
        """

        file = MOL2File(name=name)

        (geometries, all_clean_symbols, all_symbols, all_bonds, conformers) = cls._read_mol2(filename, **kwargs)
        assert len(all_bonds) == len(geometries)
        for bonds in all_bonds:
            assert isinstance(bonds, nx.Graph)
            assert len(bonds) == len(geometries[0])

        if conformers:
            # convert atom types to atomic numbers
            atomic_numbers = []
            for atom_type in all_symbols[0]:
                assert isinstance(atom_type,str), f"unexpected atom_type type: {type(atom_type)} / {atom_type}"
                fields = atom_type.split(".")
                symbol = fields[0]
                symbol = re.sub("[^A-Za-z]","",symbol)
                atomic_number = get_number(symbol)
                atomic_numbers.append(atomic_number)
            atomic_numbers = np.asarray(atomic_numbers, dtype=np.int8)

            # create ensemble
            file.ensemble = ConformationalEnsemble()
            for geometry in geometries:
                molecule = Molecule(atomic_numbers, geometry, bonds=all_bonds[0].edges, checks=False)
                file.ensemble.add_molecule(molecule, checks=False)
        else:
            file.ensemble = Ensemble()
            for this_symbols,geometry in zip(all_symbols,geometries):
                atomic_numbers=[]
                for atom_type in this_symbols:
                    assert isinstance(atom_type,str), f"unexpected atom_type type: {type(atom_type)} / {atom_type}"
                    fields = atom_type.split(".")
                    symbol = fields[0]
                    symbol = re.sub("[^A-Za-z]","",symbol)
                    atomic_number = get_number(symbol)
                    atomic_numbers.append(atomic_number)
                atomic_numbers = np.asarray(atomic_numbers, dtype=np.int8)
                molecule = Molecule(atomic_numbers, geometry, bonds=bonds.edges)
                file.ensemble.add_molecule(molecule)

        return file

    @classmethod
    def _read_mol2(
        cls, filename, contains_conformers="check", save_memory_for_conformers=True, print_status_messages=False,
    ):
        """
        Reads .mol2 files into cctk.

        Args:
            filename str): the name of the .mol2 file

            contains_conformers('check' or bool): if set to 'check', multiple geometries
                                                in the same file will be compared to see
                                                if they are conformers.  Alternatively,
                                                force the geometries to be treated as
                                                conformers (True) or not (False).  This
                                                latter option increases performance,
                                                particularly for large files.

            print_status_messages (bool): if True, update the progerss of the parsing operation to stdout.

        Returns:
            all_geometries, all_clean_symbols, all_symbols, all_bonds, contains_conformers

            all_geometries: np.ndarray(geometry number, atom number, xyz) -> position (float)
            all_clean_symbols: np.ndarray(geometry number, atom number) -> atom symbol (:obj:`str`)
            all_symbols: np.ndarray(geometry number, atom number) -> atom symbol (:obj:`str`)
            all_bonds: list(geometry_number) -> bond connectivity (:obj:`nx.Graph`)
            contains_conformers: bool (True if the geometries correspond to conformers.)
        """
        # read file
        if print_status_messages:
            print(f"Reading {filename}...", end="", flush=True)
        lines = super().read_file(filename)
        if print_status_messages:
            print(f"read {len(lines)} lines...", end="", flush=True)

        # initialize arrays
        all_geometries = []
        all_symbols = []
        all_clean_symbols = []
        all_bonds = []
        this_geometry = []
        this_symbols = []
        this_clean_symbols = []
        this_bonds = None

        # parse file
        i = 0
        in_geometry_block = False
        in_bond_block = False
        bond_number = 0
        while i < len(lines):
            # read the current line
            line = lines[i]

            # determine if we are in a geometry block
            if line.startswith("@<TRIPOS>ATOM"):
                # step forward to the first geometry line
                in_geometry_block = True
                in_bond_block = False
                i += 1
                line = lines[i]
                if contains_conformers == True and len(all_symbols) > 0:
                    this_symbols = all_symbols[0]
                    this_clean_symbols = all_clean_symbols[0]
            elif line.startswith("@<TRIPOS>BOND"):
                # update status
                in_geometry_block = False
                in_bond_block = True
                bond_number = 0

                # get next line
                i += 1
                line = lines[i]

                # initialize connectivity graph
                if len(this_geometry) == 0:
                    raise ValueError("got to bond table without a geometry")
                if contains_conformers == True and len(all_bonds) > 0:
                    this_bonds = all_bonds[0]
                else:
                    this_bonds = nx.Graph()
                    this_bonds.add_nodes_from(range(1, len(this_geometry) + 1))

            # parse geometry if appropriate
            if in_geometry_block:
                fields = line.split()
                if len(fields) < 6:
                    print("Error parsing file:")
                    print("Line = '%s'" % line.strip())
                    print(fields)
                    break
                x, y, z = float(fields[2]), float(fields[3]), float(fields[4])
                this_geometry.append([x, y, z])
                if contains_conformers != True or len(all_symbols)==0:
                    symbol = fields[5]
                    clean_symbol = fields[1]
                    this_symbols.append(symbol)
                    this_clean_symbols.append(clean_symbol)
            elif in_bond_block:
                fields = line.split()
                if len(fields) == 4 and (len(all_bonds)==0 or contains_conformers != True):
                    # parse bonds, checking that the bonds are increasing
                    try:
                        this_bond_number = int(fields[0])
                        atom1 = int(fields[1])
                        atom2 = int(fields[2])
                        n_atoms = len(this_geometry)
                        if not 1 <= atom1 <= n_atoms or not 1 <= atom2 <= n_atoms:
                            raise ValueError(f"atom number out of range: {line}")
                        if fields[3] == "ar":
                            bond_order = 1
                        else:
                            bond_order = int(fields[3])
                        if bond_order <= 0:
                            raise ValueError(f"zero or negative bond order: {line}")
                        if this_bond_number != bond_number + 1:
                            raise ValueError("non-sequential bond number")
                        bond_number = this_bond_number
                        if this_bonds.has_edge(atom1, atom2):
                            current_bond_order = this_bonds[atom1][atom2]["weight"]
                            if current_bond_order != bond_order:
                                raise ValueError(f"inconsistent bond order definition: {line}")
                        this_bonds.add_edge(atom1, atom2, weight=bond_order)
                        this_bonds.add_edge(atom2, atom1, weight=bond_order)
                    except Exception as e:
                        # assume we have left the bond block
                        in_geometry_block = False
                        in_bond_block = False
                else:
                    # we have left the bond block
                    in_geometry_block = False
                    in_bond_block = False

            # go to next line
            i += 1

            # store geometry and reinitialize if appropriate
            end_of_file = i == len(lines)
            end_of_blocks = not in_geometry_block and not in_bond_block
            if (end_of_file or end_of_blocks) and len(this_geometry) > 0:
                all_geometries.append(this_geometry)
                all_clean_symbols.append(this_clean_symbols)
                all_symbols.append(this_symbols)
                all_bonds.append(this_bonds)
                this_geometry = []
                this_symbols = []
                this_clean_symbols = []
                this_bonds = None

        # convert to numpy array
        all_geometries = np.array(all_geometries)
        all_symbols = np.array(all_symbols)
        all_clean_symbols = np.array(all_clean_symbols)

        # determine if these are conformers
        if contains_conformers == "check":
            contains_conformers = True
            for symbols, bonds in zip(all_symbols[1:], all_bonds[1:]):
                # must have the same symbols and bonds
                if not (all_symbols[0] == symbols).all() or not nx.is_isomorphic(all_bonds[0], bonds):
                    contains_conformers = False
                    break
        elif isinstance(contains_conformers, bool):
            pass
        else:
            raise ValueError("contains_conformers must be 'check' or boolean")

        # return result
        n_geometries = len(all_geometries)
        if print_status_messages:
            if n_geometries > 1:
                if contains_conformers:
                    n_atoms = len(all_geometries[0])
                    n_bonds = all_bonds[0].number_of_edges()
                    if print_status_messages:
                        print(f"read {n_geometries} conformers ({n_atoms} atoms and {n_bonds} bonds).")
                else:
                    min_n_atoms = len(all_geometries[0])
                    max_n_atoms = len(all_geometries[0])
                    for geometry in all_geometries[1:]:
                        if len(geometry) > max_n_atoms:
                            max_n_atoms = len(geometry)
                        elif len(geometry) < min_n_atoms:
                            min_n_atoms = len(geometry)
                    min_n_bonds = all_bonds[0].number_of_edges()
                    max_n_bonds = all_bonds[0].number_of_edges()
                    for bonds in all_bonds[1:]:
                        if bonds.number_of_edges() > max_n_bonds:
                            max_n_bonds = bonds.number_of_edges()
                        elif bonds.number_of_edges() < min_n_bonds:
                            min_n_bonds = bonds.number_of_edges
                    if print_status_messages:
                        print(f"read {n_geometries} unrelated geometries ({min_n_atoms}-{max_n_atoms} atoms and {min_n_bonds}-{max_n_bonds}) bonds).")
            else:
                n_atoms = len(all_geometries)
                n_bonds = all_bonds[0].number_of_edges()
                if print_status_messages:
                    print(f"read one geometry ({n_atoms} atoms and {n_bonds} bonds).")

        return (all_geometries, all_clean_symbols, all_symbols, all_bonds, contains_conformers)

    def get_molecule(self, num=None):
        """
        Returns the last molecule from the ensemble.

        If ``num`` is specified, returns ``self.ensemble.molecules[num]``
        """
        # some methods pass num=None, which overrides setting the default above
        if num is None:
            num = -1

        if not isinstance(num, int):
            raise TypeError("num must be int")

        return self.ensemble.molecules[num]

    @classmethod
    def write_molecule_to_file(cls, filename, molecule, title=None, append=False):
        """
        Write a ``.gjf`` file using the given molecule.

        Args:
            filename (str): path to the new file
            molecule (Molecule): which molecule to use -- a``Molecule`` object.
            title (str): title of the file
            append (Bool): whether to write to file normally or append
        """
        assert isinstance(molecule, Molecule), "molecule is not a valid Molecule object!"

        text = f"# {title}\n#\n#\n\n#\n#\n\n"
        text += f"@<TRIPOS>MOLECULE\n{title}\n{molecule.num_atoms()} {molecule.bonds.number_of_edges()}\nSMALL\nNO_CHARGES\n\n\n"
        text += "@<TRIPOS>ATOM\n"
        for idx, z in enumerate(molecule.atomic_numbers, start=1):
            v = molecule.get_vector(idx)
            text += f"{idx} {get_symbol(z)}{idx}    {v[0]: .4f}    {v[1]: .4f}    {v[2]: .4f} {get_symbol(z)} 0\n"
        text += "@<TRIPOS>BOND\n"
        count = 1
        for atom1, atom2, weight in molecule.bonds.edges.data("weight", default=1):
            text += f"{count} {atom1} {atom2} {weight}\n"
            count += 1

        if append:
            super().append_to_file(filename, text)
        else:
            super().write_file(filename, text)

    def write_file(self, filename, molecule=-1, **kwargs):
        """
        Write a ``.mol2`` file, using object attributes.

        Args:
            filename (str): path to the new file
            molecule (int): which molecule to use -- passed to ``self.get_molecule()``.
                Default is -1 (e.g. the last molecule), but positive integers will select from self.ensemble.molecules (0-indexed).
                A ``Molecule`` object can also be passed, in which case that molecule will be written to the file.
        """
        if molecule is None or isinstance(molecule, (np.integer, int)):
            molecule = self.ensemble.molecules[molecule]
        self.write_molecule_to_file(filename, molecule, **kwargs)

    @classmethod
    def write_ensemble_to_file(cls, filename, ensemble):
        """
        Write each structure in the specified ensemble to a single mol2 file.

        Args:
            filename (str): where to write the file
            ensemble (Ensemble): ``Ensemble`` object to write
        """
        for idx, molecule in enumerate(ensemble.molecules):
            if idx == 0:
                cls.write_molecule_to_file(filename, molecule, append=False)
            else:
                cls.write_molecule_to_file(filename, molecule, append=True)
