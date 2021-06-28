"""
Functions to handle 3D topology, graph structure, etc of ``Molecule`` objects.

Moved out of ``cctk.Molecule`` because the file was getting unwieldy.
"""

import numpy as np
import networkx as nx
import copy

from cctk.helper_functions import (
    compute_chirality,
)

def are_isomorphic(mol1, mol2, return_ordering=False):
    """
    Checks if two molecules are isomorphic (by comparing bond graphs and atomic numbers - not bond orders!).

    Args:
        mol1 (cctk.Molecule):
        mol2 (cctk.Molecule):
        return_ordering (Bool): if True, also returns a mapping between atomic numbers

    Returns:
        Boolean denoting if the molecules are isomorphic
        (optional) mapping list
    """
    assert mol1.bonds.number_of_edges() > 0, "need a bond graph to perform this operation -- try calling self.assign_connectivity()!"
    assert mol2.bonds.number_of_edges() > 0, "need a bond graph to perform this operation -- try calling self.assign_connectivity()!"

    mol1._add_atomic_numbers_to_nodes()
    mol2._add_atomic_numbers_to_nodes()

    nm = nx.algorithms.isomorphism.categorical_node_match("atomic_number", 0)
    match = nx.algorithms.isomorphism.GraphMatcher(mol1.bonds, mol2.bonds, node_match=nm)

    if match.is_isomorphic():
        if return_ordering:
            new_ordering = [match.mapping[x] for x in range(1, mol1.num_atoms() + 1)]
            return True, new_ordering
        else:
            return True
    else:
        if return_ordering:
            return False, None
        else:
            return False

def flip_meso_rings(mol, atoms):
    """
    Returns a list of permuted molecules with various ``meso`` rings renumbered.

    Args:
        mol (cctk.Molecule): molecule of interest
        atoms (list): atomic numbers of potential atoms to consider

    Returns:
        list of ``Molecule`` objects
    """
    #### get all rings in graph
    returns = [copy.deepcopy(mol)]
    for center in atoms:
        cycles = nx.cycle_basis(mol.bonds, root=center)
        for cycle in cycles:
            #### get the correct ring
            if center not in cycle:
                continue

            #### reorder to put ``center`` first
            while cycle[0] != center:
                # why yes, this /is/ a O(n) solution for reordering a list. why do you ask?
                cycle = cycle[1:] + cycle[0:1]
            assert cycle[0] == center, "graph reorder failed"

            #### create fragments
            frag1 = [cycle.pop(1)]
            frag2 = [cycle.pop(-1)]
            while len(cycle) > 2:
                frag1.append(cycle.pop(1))
                frag2.append(cycle.pop(-1))

            #### cut fragment bonds, depending on if we have even- or odd-numbered ring
            new_returns = []
            for mol in returns:
                cpy = copy.deepcopy(mol)
                cpy.remove_bond(frag1[0], cycle[0])
                cpy.remove_bond(frag2[0], cycle[0])
                if len(cycle) == 1:
                    cpy.remove_bond(frag1[-1], frag2[-1])
                elif len(cycle) == 2:
                    cpy.remove_bond(frag1[-1], cycle[-1])
                    cpy.remove_bond(frag2[-1], cycle[-1])

                #### generate graphs
                graph1 = None
                graph2 = None
                fragments = nx.connected_components(cpy.bonds)
                for fragment in fragments:
                    if frag1[0] in fragment:
                        graph1 = cpy.bonds.subgraph(fragment)
                    if frag2[0] in fragment:
                        graph2 = cpy.bonds.subgraph(fragment)

                assert isinstance(graph1, nx.Graph), "can't find graph 1"
                assert isinstance(graph2, nx.Graph), "can't find graph 1"

                #### do our two ring-halves match?? if so, we swap them
                nm = nx.algorithms.isomorphism.categorical_node_match("atomic_number", 0)
                match = nx.algorithms.isomorphism.GraphMatcher(graph1, graph2, node_match=nm)

                if match.is_isomorphic():
                    for k,v in match.mapping.items():
                        cpy = cpy.swap_atom_numbers(k, v)

                    #### redo all the bonds we ablated
                    if len(cycle) == 1:
                        cpy.add_bond(frag1[-1], frag2[-1], mol.get_bond_order(frag1[-1], frag2[-1]))
                    elif len(cycle) == 2:
                        cpy.add_bond(frag1[-1], cycle[-1], mol.get_bond_order(frag1[-1], cycle[-1]))
                        cpy.add_bond(frag2[-1], cycle[-1], mol.get_bond_order(frag2[-1], cycle[-1]))
                    cpy.add_bond(frag1[0], cycle[0], mol.get_bond_order(frag1[0], cycle[0]))
                    cpy.add_bond(frag2[0], cycle[0], mol.get_bond_order(frag2[0], cycle[0]))

                    new_returns.append(cpy)
            returns = returns + new_returns
    return returns

def exchange_identical_substituents(mol, center, self_permutations=None):
    """
    Replace homotopic/enantiotopic/diastereotopic substituents about a single atom.

    If a list of permuted ``Molecule`` objects is passed (as ``self_permutations``), then this code will apply this to each member and return a list.

    Args:
        mol (cctk.Molecule): molecule of interest
        center (integer): atomic number of atom to swap substituents around
        self_permutations (list of Molecules): optional list of starting ``Molecule`` objects

    Returns:
        ``Molecule`` object (or list if ``self_permutations`` is not ``None``)
    """
    assert mol.bonds.number_of_edges() > 0, "need a bond graph to perform this operation -- try calling self.assign_connectivity()!"
    mol._add_atomic_numbers_to_nodes()
    neighbors = list(mol.bonds[center])

    returns = [copy.deepcopy(mol)]
    if self_permutations is not None:
        returns = self_permutations


    for i in range(len(neighbors)):
        for j in range(i+1, len(neighbors)):
            try:
                _, frag1 = mol._get_bond_fragments(center, neighbors[i])
                _, frag2 = mol._get_bond_fragments(center, neighbors[j])

                graph1 = mol.bonds.subgraph(frag1)
                graph2 = mol.bonds.subgraph(frag2)

                nm = nx.algorithms.isomorphism.categorical_node_match("atomic_number", 0)
                match = nx.algorithms.isomorphism.GraphMatcher(graph1, graph2, node_match=nm)
                if match.is_isomorphic():
                    for m in returns:
                        new_mol = copy.deepcopy(m)
                        for k,v in match.mapping.items():
                            new_mol = new_mol.swap_atom_numbers(k, v)
                        if self_permutations is None:
                            return new_mol

                    returns.append(new_mol)

            except ValueError as e:
                pass # probably indicates a cycle

    if self_permutations is None:
        raise ValueError("could not find substituents to switch")
    else:
        return returns

def get_chirality_report(mol, centers=None):
    """
    Computes chirality at stereogenic centers.

    Args:
        mol (cctk.Molecule): molecule of interest
        centers (list): atomic numbers to check. defaults to all centers with 4+ substituents.

    Returns:
        dict with centers as keys and ±1 as values
    """
    if centers is None:
        centers = get_stereogenic_centers(mol)
    assert isinstance(centers, list)

    results = {}
    for center in centers:
        neighbors = list(mol.bonds[center])
        neighbors.sort()
        assert len(neighbors) >= 4, f"atom {center} has fewer than 4 neighbors ({neighbors})!"
        results[center] = compute_chirality(*[mol.get_vector(n, center) for n in neighbors])

    return results

def get_stereogenic_centers(mol):
    """
    Returns every atom making 4 or more bonds. A bit misleading, since diastereotopic protons/meso protons are also counted.
    """
    assert mol.bonds.number_of_edges() > 0, "need a bond graph to perform this operation -- try calling self.assign_connectivity()!"
    num_neighbors = np.array([len(list(mol.bonds[x])) for x in range(1, mol.num_atoms() + 1)])
    return [int(x) for x in list(np.ravel(np.argwhere(num_neighbors >= 4)) + 1)] # love me some off-by-one indexing errors

def get_exchangeable_centers(mol):
    """
    Returns all atoms making 4 or more bonds that have two isomorphic substituents, i.e. where renumbering could be broken.
    """
    centers = get_stereogenic_centers(mol)
    exchangeable_centers = []
    for center in centers:
        try:
            exchange_identical_substituents(mol, center)
            exchangeable_centers.append(center)
            continue
        except Exception as e:
            pass

        mols = flip_meso_rings(mol, atoms=[center])
        if len(mols) > 1:
            exchangeable_centers.append(center)

    return exchangeable_centers

def find_group(mol, group):
    """
    Finds instances of ``group`` within ``mol``.

    Args:
        mol (cctk.Molecule): molecule to search within
        group (cctk.Group): group to search for

    Returns:
        list of dictionaries mapping from molecule atomic numbers to group atomic numbers
    """
    assert mol.bonds.number_of_edges() > 0, "need a bond graph to perform this operation -- try calling self.assign_connectivity()!"
    assert group.bonds.number_of_edges() > 0, "need a bond graph to perform this operation -- try calling self.assign_connectivity()!"

    mol._add_atomic_numbers_to_nodes()
    group._add_atomic_numbers_to_nodes()
    group_map = group.map_from_truncated()
    group.remove_atom(group.attach_to)

    nm = nx.algorithms.isomorphism.categorical_node_match("atomic_number", 0)
    match = nx.algorithms.isomorphism.GraphMatcher(mol.bonds, group.bonds, node_match=nm)

    #### need to only find unique mappings - combinations, not permutations
    mappings = []
    for sg in match.subgraph_isomorphisms_iter():
        unique = True
        for m in mappings:
            if set(m.keys()) == set(sg.keys()):
                unique = False
                break
        if unique:
            mappings.append(sg)

    composition = [{k: group_map[v] for k, v in m.items()} for m in mappings]
    return composition

