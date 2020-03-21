import re
import numpy as np
import networkx as nx

from abc import abstractmethod

from cctk import File, Ensemble, ConformationalEnsemble
from cctk.helper_functions import get_symbol, get_number


class MOL2File(File):
    """
    Generic class for all ``.mol2`` files.

    Attributes:
        name (str): name of file
        molecules (Ensemble): ``Ensemble`` or ``ConformationalEnsemble`` object
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

        (geometries, symbols, atom_types, bonds, conformers) = cls._read_mol2(filename, **kwargs)
        atomic_numbers = np.array([get_number(z) for z in symbols], dtype=np.int8)
        atomic_numbers = np.tile(atomic_numbers, (len(geometries), 1))
        if conformers == True:
            file.molecules = ConformationalEnsemble(geometries=geometries, atomic_numbers=atomic_numbers, bonds=[bonds.edges() for g in geometries])
        else:
            file.molecules = Ensemble(geometries=geometries, atomic_numbers=atomic_numbers, bonds=[b.edges() for b in bonds])

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

            save_memory_for_conformers (bool): if True, the first dimension (geometry number) will be
                                            dropped from the symbols and bonds to prevent
                                            the storage of redundant information.  Thus,
                                            symbols will be a one-dimensional :obj:`np.array` of
                                            :obj:`str` and bonds will be a single :obj:`nx.Graph`.

            print_status_messages (bool): if True, update the progerss of the parsing operation to stdout.

        Returns:
            all_geometries, all_clean_symbols, all_symbols, all_bonds, contains_conformers

            all_geometries: np.array(geometry number, atom number, xyz) -> position (float)
            all_clean_symbols: np.array(geometry number, atom number) -> atom symbol (:obj:`str`)
            all_symbols: np.array(geometry number, atom number) -> atom symbol (:obj:`str`)
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
                this_bonds = nx.Graph()
                this_bonds.add_nodes_from(range(1, len(this_geometry) + 1))

            # parse geometry if appropriate
            if in_geometry_block:
                fields = re.split(" +", line.strip())
                if len(fields) < 6:
                    print("Error parsing file:")
                    print("Line = '%s'" % line.strip())
                    print(fields)
                    break
                x, y, z = float(fields[2]), float(fields[3]), float(fields[4])
                this_geometry.append([x, y, z])
                symbol = fields[5]
                clean_symbol = fields[1]
                this_symbols.append(symbol)
                this_clean_symbols.append(clean_symbol)
            elif in_bond_block:
                fields = re.split(" +", line.strip())
                if len(fields) == 4:
                    # parse bonds, checking that the bonds are increasing
                    try:
                        this_bond_number = int(fields[0])
                        atom1 = int(fields[1])
                        atom2 = int(fields[2])
                        n_atoms = len(this_geometry)
                        if not 1 <= atom1 <= n_atoms or not 1 <= atom2 <= n_atoms:
                            raise ValueError(f"atom number out of range: {line}")
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
                    except:
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

        # if requested, just store one copy of all_symbols and all_bonds
        if save_memory_for_conformers and contains_conformers:
            all_symbols = all_symbols[0]
            all_clean_symbols = all_clean_symbols[0]
            all_bonds = all_bonds[0]

        # return result
        n_geometries = len(all_geometries)
        if print_status_messages:
            if n_geometries > 1:
                if contains_conformers:
                    n_atoms = len(all_geometries[0])
                    n_bonds = all_bonds.number_of_edges()
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
                n_bonds = all_bonds.number_of_edges()
                if print_status_messages:
                    print(f"read one geometry ({n_atoms} atoms and {n_bonds} bonds).")

        #### sometimes these labels switch? so gotta check
        if len(all_clean_symbols[0]) < len(all_symbols[0]):
            return (all_geometries, all_clean_symbols, all_symbols, all_bonds, contains_conformers)
        else:
            return (all_geometries, all_symbols, all_clean_symbols, all_bonds, contains_conformers)
