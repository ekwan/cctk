import unittest, sys, os, io, copy
import numpy as np
import cctk

import cctk.helper_functions as helper

class TestElement(unittest.TestCase):
    def test_get_number(self):
        self.assertEqual(helper.get_number("CL"), 17)
        self.assertEqual(helper.get_number("cl"), 17)
        self.assertEqual(helper.get_number("Cl"), 17)

        self.assertEqual(helper.get_number("Bq"), 0)
        self.assertEqual(helper.get_number("U"), 92)

    def test_get_symbol(self):
        self.assertEqual(helper.get_symbol(0), "Bq")
        self.assertEqual(helper.get_symbol(62), "Sm")
        self.assertEqual(helper.get_symbol(46), "Pd")
        self.assertEqual(helper.get_symbol(6), "C")

if __name__ == '__main__':
    unittest.main()
