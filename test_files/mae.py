import re
import numpy as np
import networkx as nx

### Mae File Format Reader ###
# reads uncompressed macromodel files
#
# returns:
# geometries: np.array(geometry number, atom number, xyz) -> position (float)
# symbols:    np.array(geometry number, atom number) -> symbol (str)
# bonds:      np.array(geometry number) -> connectivity (nx.Graph)
# property_names: np.array(geometry number) -> str
# property_values: np.array(conformer, property) -> property_value (str)
# contains_conformers: whether this file contains conformers (bool)
def read_mae(
    filename,
    contains_conformers="check",
    save_memory_for_conformers=True,
    print_status_messages=True,
):
    # read file
    if print_status_messages:
        print(f"Reading {filename}...", end="", flush=True)
    with open(filename, "r") as filehandle:
        lines = filehandle.read().splitlines()
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
                    raise ValueError(
                        f"expected ':::' here but line {i+1} is:\n{next_line}\n"
                    )
                current_block_type = "property_values"
                i += 1
            elif line.startswith(":::"):
                raise ValueError(
                    f"expected to see i_m_ct_format as the last property (line {i+1})"
                )
            else:
                fields = re.split(" +", line)
                if len(fields) != 1:
                    raise ValueError(
                        f"unexpected number of fields in property name line: {line}"
                    )
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
                    raise ValueError(
                        f"non-sequential bond number (expected {this_bonds.number_of_edges()+1} but got {bond_number})"
                    )
                if this_bonds.has_edge(atom1, atom2):
                    current_bond_order = this_bonds[atom1][atom2]["weight"]
                    if current_bond_order != bond_order:
                        raise ValueError(f"inconsistent bond order definition: {line}")
                this_bonds.add_edge(atom1, atom2, weight=bond_order)
                this_bonds.add_edge(atom2, atom1, weight=bond_order)
        # print(f"{i+1:6d} : {current_block_type} : {line}")

    # return result
    if print_status_messages:
        print("done.")
    return (
        geometries,
        symbols,
        bonds,
        property_names,
        property_values,
        contains_conformers,
    )
