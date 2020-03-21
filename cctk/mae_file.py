import re
import numpy as np
import networkx as nx

from abc import abstractmethod

from cctk import File, Ensemble, ConformationalEnsemble
from cctk.helper_functions import get_symbol, get_number


class MAEFile(File):
    """
    Generic class for all ``.mae`` files.

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
        Reads ``.mae`` file and generates a ``MAEFile`` instance.

        Args:
            filename (str): path to file
            name (str): name of the file

        Returns:
            MAEFile object
            property names (list)
            property_values (list)
        """

        file = MAEFile(name=name)

        (geometries, symbols, bonds, p_names, p_vals, conformers) = cls._read_mae(filename, **kwargs)
        atomic_numbers = np.array([get_number(z) for z in symbols], dtype=np.int8)
        atomic_numbers = np.tile(atomic_numbers, (len(geometries), 1))
        if conformers == True:
            file.molecules = ConformationalEnsemble(geometries=geometries, atomic_numbers=atomic_numbers, bonds=[bonds.edges() for g in geometries])
        else:
            file.molecules = Ensemble(geometries=geometries, atomic_numbers=atomic_numbers, bonds=[b.edges() for b in bonds])

        return file, p_names, p_vals

    @classmethod
    def _read_mae(
        cls, filename, contains_conformers="check", save_memory_for_conformers=True, print_status_messages=False,
    ):
        """
        Reads uncompressed Macromodel files.

        Args:
            filename (str): path to file
            contains_conformers (str): one of ``check``, ``True``, or ``False``
            save_memory_for_conformers (Bool):
            print_status_messages (Bool):

        Returns:
            geometries (np.array): array of 3-tuples of geometries
            symbols (np.array): array of atom symbols (str)
            bonds (nx.Graph): ``NetworkX`` graph of bond information
            property_names:
            property_values:
            contains_conformers (Bool): whether or not the file contains conformers
        """
        # read file
        if print_status_messages:
            print(f"Reading {filename}...", end="", flush=True)
        lines = super().read_file(filename)
        if print_status_messages:
            print(f"read {len(lines)} lines...", end="", flush=True)

        # initialize arrays
        geometries = []
        symbols = []
        bonds = []
        property_names = []
        property_values = []
        this_geometry = None
        this_symbols = None
        this_bonds = None
        this_property_names = None
        this_property_values = None

        # parse file
        i = 0
        current_block_type = None
        while i < len(lines):
            # read the current line
            line = lines[i]
            i += 1

            # determine if we are in a molecule block
            end_of_file = i + 1 == len(lines)
            if current_block_type is None and (line.startswith("f_m_ct") or end_of_file):
                # store the current results if any
                if this_geometry is not None and len(this_geometry) > 0:
                    geometries.append(this_geometry)
                    symbols.append(this_symbols)
                    bonds.append(this_bonds)
                    property_names.append(this_property_names)
                    property_values.append(this_property_values)

                # prepare to read a new molecule
                current_block_type = "property_names"
                this_geometry = []
                this_symbols = []
                this_bonds = None
                this_property_names = []
                this_property_values = []
                continue

            # read property names
            elif current_block_type == "property_names":
                line = line.strip()
                if line.startswith("i_m_ct_format"):
                    next_line = lines[i].strip()
                    if next_line != ":::":
                        raise ValueError(f"expected ':::' here but line {i+1} is:\n{next_line}\n")
                    current_block_type = "property_values"
                    i += 1
                elif line.startswith(":::"):
                    raise ValueError(f"expected to see i_m_ct_format as the last property (line {i+1})")
                else:
                    fields = re.split(" +", line)
                    if len(fields) != 1:
                        raise ValueError(f"unexpected number of fields in property name line: {line}")
                    this_property_names.append(line)

            # read property values
            elif current_block_type == "property_values":
                n_properties = len(this_property_names)
                for j in range(n_properties):
                    this_property_values.append(lines[i + j])
                i += n_properties
                current_block_type = "looking_for_geometry1"

            # look for geometry block
            elif current_block_type == "looking_for_geometry1":
                if line.startswith("  m_atom"):
                    current_block_type = "looking_for_geometry2"
            elif current_block_type == "looking_for_geometry2":
                if line.strip() == ":::":
                    current_block_type = "geometry_block"

            # parse geometry
            elif current_block_type == "geometry_block":
                line = line.strip()
                if line == ":::":
                    current_block_type = "bond_block"

                    # initialize bond connectivity graph
                    this_bonds = nx.Graph()
                    n_atoms = len(this_symbols)
                    this_bonds.add_nodes_from(range(1, n_atoms + 1))
                    i += 7
                else:
                    fields = re.split(" +", line)
                    x, y, z = float(fields[2]), float(fields[3]), float(fields[4])
                    this_geometry.append((x, y, z))
                    symbol = fields[-1]
                    this_symbols.append(symbol)

            # parse bonds
            elif current_block_type == "bond_block":
                line = line.strip()
                if line == ":::":
                    current_block_type = None
                else:
                    fields = re.split(" +", line)
                    bond_number, atom1, atom2, bond_order = (
                        int(fields[0]),
                        int(fields[1]),
                        int(fields[2]),
                        int(fields[3]),
                    )
                    n_atoms = len(this_geometry)
                    if not 1 <= atom1 <= n_atoms or not 1 <= atom2 <= n_atoms:
                        raise ValueError(f"atom number out of range: {line}")
                    bond_order = int(fields[3])
                    if bond_order <= 0:
                        raise ValueError(f"zero or negative bond order: {line}")
                    if this_bonds.number_of_edges() != bond_number - 1:
                        raise ValueError(f"non-sequential bond number (expected {this_bonds.number_of_edges()+1} but got {bond_number})")
                    if this_bonds.has_edge(atom1, atom2):
                        current_bond_order = this_bonds[atom1][atom2]["weight"]
                        if current_bond_order != bond_order:
                            raise ValueError(f"inconsistent bond order definition: {line}")
                    this_bonds.add_edge(atom1, atom2, weight=bond_order)
                    this_bonds.add_edge(atom2, atom1, weight=bond_order)

        # convert to numpy array
        geometries = np.array(geometries)
        symbols = np.array(symbols)
        property_names = np.array(property_names)
        property_values = np.array(property_values)

        # determine if these are conformers
        if contains_conformers == "check":
            contains_conformers = True
            for this_symbols, this_bonds in zip(symbols[1:], bonds[1:]):
                # must have the same symbols and bonds
                if not (symbols[0] == this_symbols).all() or not nx.is_isomorphic(bonds[0], this_bonds):
                    contains_conformers = False
                    break
        elif isinstance(contains_conformers, bool):
            pass
        else:
            raise ValueError("contains_conformers must be 'check' or boolean")

        # if requested, just store one copy of symbols and bonds
        if save_memory_for_conformers and contains_conformers:
            symbols = symbols[0]
            bonds = bonds[0]

        # return result
        n_geometries = len(geometries)
        if print_status_messages:
            if n_geometries > 1:
                if contains_conformers:
                    n_atoms = len(geometries[0])
                    n_bonds = bonds.number_of_edges()
                    if print_status_messages:
                        print(f"read {n_geometries} conformers ({n_atoms} atoms and {n_bonds} bonds).")
                else:
                    min_n_atoms = len(geometries[0])
                    max_n_atoms = len(geometries[0])
                    for geometry in geometries[1:]:
                        if len(geometry) > max_n_atoms:
                            max_n_atoms = len(geometry)
                        elif len(geometry) < min_n_atoms:
                            min_n_atoms = len(geometry)
                    min_n_bonds = bonds[0].number_of_edges()
                    max_n_bonds = bonds[0].number_of_edges()
                    for this_bonds in bonds[1:]:
                        if this_bonds.number_of_edges() > max_n_bonds:
                            max_n_bonds = this_bonds.number_of_edges()
                        elif this_bonds.number_of_edges() < min_n_bonds:
                            min_n_bonds = bonds.number_of_edges
                    if print_status_messages:
                        print(f"read {n_geometries} unrelated geometries ({min_n_atoms}-{max_n_atoms} atoms and {min_n_bonds}-{max_n_bonds}) bonds).")
            else:
                n_atoms = len(geometries)
                n_bonds = bonds.number_of_edges()
                if print_status_messages:
                    print(f"read one geometry ({n_atoms} atoms and {n_bonds} bonds).")

        # return result
        return (
            geometries,
            symbols,
            bonds,
            property_names,
            property_values,
            contains_conformers,
        )
