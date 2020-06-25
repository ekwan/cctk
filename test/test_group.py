import unittest, sys, os, io, copy
import numpy as np
import cctk
import cctk.topology as top

from cctk.load_groups import load_group

class TestGroup(unittest.TestCase):
    def test_group_add(self):
        path = "test/static/acetaldehyde.out"
        old_path = "test/static/14-butanedione.gjf"
        new_path = "test/static/new_14-butanedione.gjf"

        file = cctk.GaussianFile.read_file(path)
        group = cctk.Group.new_from_molecule(attach_to=6, molecule=file.get_molecule())
        new_mol, mol_map, group_map = cctk.Group.add_group_to_molecule(file.get_molecule(), group, 5, return_mapping=True)
        file.write_file("test/static/new_14-butanedione.gjf", molecule=new_mol)

        # test mapping
        for z in range(1, file.get_molecule().num_atoms()+1):
            if mol_map[z] is not None:
                self.assertEqual(file.get_molecule().atomic_numbers[z], new_mol.atomic_numbers[mol_map[z]])

        for z in range(1, group.num_atoms()+1):
            if group_map[z] is not None:
                self.assertEqual(group.atomic_numbers[z], new_mol.atomic_numbers[group_map[z]])

        # test geometry
        with open(old_path) as old:
            with open(new_path) as new:
                self.assertListEqual(list(old), list(new))

        os.remove(new_path)

    def test_group_remove(self):
        path = "test/static/14-butanedione.gjf"

        file = cctk.GaussianFile.read_file(path)
        mol = file.get_molecule().assign_connectivity()

        mol2, group, map_m, map_g = cctk.Group.remove_group_from_molecule(mol, 2, 5, return_mapping=True)

        self.assertTrue(isinstance(mol2, cctk.Molecule))
        self.assertTrue(isinstance(group, cctk.Group))

        for x in [mol2, group]:
            self.assertEqual(x.num_atoms(), 7)

        for z in range(1, mol.num_atoms()+1):
            if z == 2 or z == 5:
                continue
            if z in map_m.keys() and map_m[z] is not None:
                self.assertEqual(mol.atomic_numbers[z], mol2.atomic_numbers[map_m[z]])
            if z in map_g.keys() and map_g[z] is not None:
                self.assertEqual(mol.atomic_numbers[z], group.atomic_numbers[map_g[z]])

    def test_epimerize(self):
        path = "test/static/L-Ala.gjf"
        new_path = "test/static/L-Ala2.gjf"

        file = cctk.GaussianFile.read_file(path)
        l_ala = file.get_molecule().assign_connectivity()

        d_ala = l_ala.epimerize(1, 3, 11)
        self.assertListEqual(list(l_ala.atomic_numbers), list(d_ala.atomic_numbers))

        l_ala2 = d_ala.epimerize(1, 3, 11)
        file.write_file(new_path, l_ala2)

        for x1, x2 in zip(np.ravel(l_ala.geometry), np.ravel(l_ala2.geometry)):
            self.assertTrue(abs(float(x1)-float(x2)) < 0.0001)

        os.remove(new_path)

    def test_subgroup_isomorphism(self):
        group = load_group("methyl")
        molecule = cctk.Molecule.new_from_name("camphor").assign_connectivity()

        subgroups = top.find_group(molecule, group)
        self.assertEqual(len(subgroups), 3)

        molecule = cctk.Molecule.new_from_name("tert-butanol").assign_connectivity()
        sets = molecule.get_symmetric_atoms()
        self.assertEqual(len(sets[0]), 3)
        self.assertEqual(len(sets[1]), 9)

if __name__ == '__main__':
    unittest.main()
