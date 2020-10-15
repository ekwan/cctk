import unittest, sys, os, io, copy, math
import numpy as np
import random
import networkx as nx

import cctk
import cctk.topology as top

class TestRenumber(unittest.TestCase):

    def test_renumber(self):
        # load test files
        original_molecules = {}
        n_original_molecules = 8
        for i in range(n_original_molecules):
            filename = f"test/static/renumber_{i}.gjf"
            gaussian_file = cctk.GaussianFile.read_file(filename)
            molecule = gaussian_file.get_molecule()
            molecule.assign_connectivity(cutoff=0.1)
            original_molecules[i] = molecule

        # function for checking if two molecules have the same numbering
        def check_numbering(m1, m2, comparison_name="comparison"):
            n1 = m1.get_n_atoms()
            n2 = m2.get_n_atoms()
            self.assertTrue(n1 == n2, f"molecules have a different number of atoms: {n1} vs. {n2}")

            elements1 = m1.atomic_numbers
            elements2 = m2.atomic_numbers
            for i,(e1,e2) in enumerate(zip(elements1, elements2)):
                self.assertTrue(e1 == e2, f"atomic numbers do not match for atom {i+1}: {e1} vs. {e2}")

            g1 = m1.bonds
            g2 = m2.bonds
            for atom_number in range(1,n1+1):
                neighbors1 = list(sorted(g1.neighbors(atom_number)))
                neighbors2 = list(sorted(g2.neighbors(atom_number)))
                self.assertListEqual(neighbors1, neighbors2, f"\ncomparison {comparison_name}: neighbors for atom {atom_number} don't match: {neighbors1} vs. {neighbors2}")

        # function for permuting a molecule
        def permute(molecule, n_shuffles=50):
            molecule2 = copy.deepcopy(molecule)
            n_atoms = molecule2.get_n_atoms()
            molecule2.bonds = nx.Graph()
            molecule2.bonds.add_nodes_from(range(1, n_atoms + 1))
            choices = list(range(1, n_atoms))
            atomic_numbers = molecule2.atomic_numbers
            geometry = molecule2.geometry
            for i in range(n_shuffles):
                a1,a2 = random.sample(choices, 2)
                atomic_numbers[a1], atomic_numbers[a2] = atomic_numbers[a2], atomic_numbers[a1]
                positions_a1, positions_a2 = geometry[a1].copy(), geometry[a2].copy()
                geometry[a1] = positions_a2
                geometry[a2] = positions_a1
            molecule2.assign_connectivity(cutoff=0.1)
            return molecule2

        # check original molecules
        for i in range(n_original_molecules-1):
            for j in range(i,n_original_molecules):
                # also checks a molecule against itself for extra safety
                m1 = original_molecules[i]
                m2 = original_molecules[j]
                check_numbering(m1, m2, f"orig {i} vs. {j}")

        # permute and check renumbering
        template_molecule = original_molecules[0]
        n_permutations_per_molecule = 1
        for i in range(1,n_original_molecules):
            original_molecule = original_molecules[i]
            for j in range(n_permutations_per_molecule):
                permuted_molecule = permute(original_molecule)
                try:
                    renumbered_molecule = permuted_molecule.renumber_to_match(template_molecule)
#                    renumbered_molecule = permuted_molecule.renumber_to_match(original_molecule)
                    cctk.GaussianFile.write_molecule_to_file("test/static/renumber_error_4_renumbered.gjf", renumbered_molecule, "#p")
                    check_numbering(original_molecule, renumbered_molecule, f"mol. {i} vs. permuted")
                except Exception as e:
                    cctk.GaussianFile.write_molecule_to_file("test/static/renumber_error_1_template.gjf", template_molecule, "#p")
                    cctk.GaussianFile.write_molecule_to_file("test/static/renumber_error_2_original.gjf", original_molecule, "#p")
                    cctk.GaussianFile.write_molecule_to_file("test/static/renumber_error_3_permuted.gjf", permuted_molecule, "#p")
#                    cctk.GaussianFile.write_molecule_to_file("test/static/renumber_error_4_renumbered.gjf", renumbered_molecule, "#p")
                    self.assertTrue(False, f"error trying to renumber:\n{e}")

    def test_diastereotopicity(self):
        mol1 = cctk.GaussianFile.read_file("test/static/chiral_fluorine.gjf").get_molecule().assign_connectivity()
        mol2 = cctk.GaussianFile.read_file("test/static/chiral_fluorine_2.gjf").get_molecule().assign_connectivity() # two protons switched (10 & 11)

        self.assertTrue(np.array_equal(mol1.atomic_numbers, mol2.atomic_numbers))
        self.assertFalse(np.array_equal(mol1.geometry, mol2.geometry))

        mol2 = mol2.renumber_to_match(mol1)
        self.assertTrue(np.array_equal(mol1.atomic_numbers, mol2.atomic_numbers))
        self.assertTrue(np.array_equal(mol1.geometry, mol2.geometry))

        mol3 = cctk.GaussianFile.read_file("test/static/chiral_fluorine_3.gjf").get_molecule().assign_connectivity() # just a bond rotation
        self.assertListEqual(list(top.get_chirality_report(mol1).values()), list(top.get_chirality_report(mol3).values()))

        mol4 = cctk.GaussianFile.read_file("test/static/chiral_fluorine_4.gjf").get_molecule().assign_connectivity() # two protons switched (24 & 25)
        self.assertTrue(np.array_equal(mol1.atomic_numbers, mol4.atomic_numbers))
        self.assertFalse(np.array_equal(mol1.geometry, mol4.geometry))

        mol4 = mol4.renumber_to_match(mol1)
        self.assertTrue(np.array_equal(mol1.atomic_numbers, mol4.atomic_numbers))
        self.assertTrue(np.array_equal(mol1.geometry, mol4.geometry))

        mol5 = cctk.GaussianFile.read_file("test/static/chiral_fluorine_5.gjf").get_molecule().assign_connectivity() # two carbons switched (16 & 17)
        self.assertTrue(np.array_equal(mol1.atomic_numbers, mol5.atomic_numbers))
        self.assertFalse(np.array_equal(mol1.geometry, mol5.geometry))

        mol5 = mol5.renumber_to_match(mol1)
        self.assertTrue(np.array_equal(mol1.atomic_numbers, mol5.atomic_numbers))
        self.assertTrue(np.array_equal(mol1.geometry, mol5.geometry))

        mol6 = cctk.GaussianFile.read_file("test/static/chiral_fluorine_6.gjf").get_molecule().assign_connectivity()
        mol7 = cctk.GaussianFile.read_file("test/static/chiral_fluorine_7.gjf").get_molecule().assign_connectivity() # switched gem-dimethyls
        self.assertTrue(np.array_equal(mol6.atomic_numbers, mol7.atomic_numbers))
        self.assertFalse(np.array_equal(mol6.geometry, mol7.geometry))

        mol6 = mol6.renumber_to_match(mol7)
        self.assertTrue(np.array_equal(mol6.atomic_numbers, mol7.atomic_numbers))
        self.assertTrue(np.array_equal(mol6.geometry, mol7.geometry))

        mol8 = cctk.GaussianFile.read_file("test/static/chiral_fluorine_8.gjf").get_molecule().assign_connectivity() # meso cyclohexane ring added
        mol9 = cctk.GaussianFile.read_file("test/static/chiral_fluorine_9.gjf").get_molecule().assign_connectivity() # reorgnized
        self.assertTrue(np.array_equal(mol8.atomic_numbers, mol9.atomic_numbers))
        self.assertFalse(np.array_equal(mol8.geometry, mol9.geometry))

        mol8 = mol8.renumber_to_match(mol9)
        self.assertTrue(np.array_equal(mol8.atomic_numbers, mol9.atomic_numbers))
        self.assertTrue(np.array_equal(mol8.geometry, mol9.geometry))

if __name__ == '__main__':
    unittest.main()
