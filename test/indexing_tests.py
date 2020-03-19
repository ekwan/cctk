import unittest, sys, os, io, copy
import numpy as np
import cctk

class TestOneIndexedArray(unittest.TestCase):
    def test_indexing(self):
        array = [1,2,3,4,5]
        new_a = cctk.OneIndexedArray(array)

        self.assertEqual(new_a[5], 5)
        self.assertEqual(new_a[1], 1)
        new_a[1] = 8
        self.assertEqual(new_a[1], 8)

        a_2d = [[1,1,1],[2,2,2],[3,3,3]]
        new_a_2d = cctk.OneIndexedArray(a_2d)
        self.assertEqual(new_a_2d[1,1], 1)
        self.assertEqual(new_a_2d[2,0], 2)

        new_a_2d[3,0] = 7
        self.assertEqual(new_a_2d[3,0], 7)

if __name__ == '__main__':
    unittest.main()
