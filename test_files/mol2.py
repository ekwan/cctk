import re
import numpy as np
import networkx as nx

### Mol2 File Format Reader ###
# reads a mol2 file
#
# returns:
# geometries: np.array(conformer, atom number, xyz)
# symbols:  np.array(conformer, atom number) -> symbol (str)
# bonds:    np.array(conformer) -> nx.graph
#
# If save_memory_for_conformers is True, the conformer dimension will be dropped
# from geometries, symbols, and bonds.  Additionally, only one copy of the symbols
# and bonds will be stored.  That is, symbols will be just an np.array of strings
# and bonds will be just a single nx.Graph.
#
# contains_conformers can be 'check', True, or False
def read_mol2(filename, contains_conformers='check', save_memory_for_conformers=True, print_status_messages=True):
    # read file
    if print_status_messages:
        print(f"Reading {filename}...", end="", flush=True)
    with open(filename, "r") as filehandle:
        lines = filehandle.read().splitlines()
    if print_status_messages:
        print(f"read {len(lines)} lines...", end="", flush=True)

    # initialize arrays
    all_geometries = []
    all_symbols = []
    all_bonds = []
    this_geometry = []
    this_symbols = []
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
            this_symbols.append(symbol)
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
                            raise ValueError(
                                f"inconsistent bond order definition: {line}"
                            )
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
        end_of_blocks = (not in_geometry_block and not in_bond_block)
        if ( end_of_file or end_of_blocks ) and len(this_geometry) > 0:
            all_geometries.append(this_geometry)
            all_symbols.append(this_symbols)
            all_bonds.append(this_bonds)
            this_geometry = []
            this_symbols = []
            this_bonds = None

    # convert to numpy array
    all_geometries = np.array(all_geometries)
    all_symbols = np.array(all_symbols)

    # determine if these are conformers
    if contains_conformers == 'check':
        contains_conformers = True
        for symbols, bonds in zip(all_symbols[1:], all_bonds[1:]):
            # must have the same symbols and bonds
            if not (all_symbols[0] == symbols).all() or not nx.is_isomorphic(all_bonds[0], bonds):
                contains_conformers = False
                break
    elif isinstance(contains_conformers,bool):
        pass
    else:
        raise ValueError("contains_conformers must be 'check' or boolean")

    # if requested, just store one copy of all_symbols and all_bonds
    if save_memory_for_conformers and contains_conformers:
        all_symbols = all_symbols[0]
        all_bonds = all_bonds[0]

    # return result
    n_geometries = len(all_geometries)
    if print_status_messages:
        if n_geometries > 1:
            if contains_conformers:
                n_atoms = len(all_geometries[0])
                n_bonds = all_bonds.number_of_edges()
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
                print(
                    f"read {n_geometries} unrelated geometries ({min_n_atoms}-{max_n_atoms} atoms and {min_n_bonds}-{max_n_bonds}) bonds)."
                )
        else:
            n_atoms = len(all_geometries)
            n_bonds = all_bonds.number_of_edges()
            print(f"read one geometry ({n_atoms} atoms and {n_bonds} bonds).")
    return (all_geometries, all_symbols, all_bonds, contains_conformers)

