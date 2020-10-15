import unittest, sys, os, io, copy
import numpy as np
import cctk

#### tests indexing for the OneIndexedArray object
class TestOneIndexedArray(unittest.TestCase):
    def test_indexing(self):
        array = [1,2,3,4,5]
        new_a = cctk.OneIndexedArray(array)

        self.assertEqual(len(new_a), 5)
        self.assertEqual(str(new_a), "[1 2 3 4 5]")
        self.assertEqual(str(new_a[1:3]), "[1 2]")
        self.assertEqual(new_a[5], 5)
        self.assertEqual(new_a[1], 1)

        self.assertTrue(isinstance(cctk.OneIndexedArray(new_a), cctk.OneIndexedArray))

        self.assertEqual(new_a[[1]], 1)
        self.assertListEqual(list(new_a[[1,2]]), [1,2])

        new_a[1] = 8
        self.assertEqual(str(new_a), "[8 2 3 4 5]")
        self.assertEqual(new_a[1], 8)

        a_2d = [[1,1,1],[2,2,2],[3,3,3]]

        new_a_2d = cctk.OneIndexedArray(a_2d)
        self.assertEqual(new_a_2d[1,1], 1)
        self.assertEqual(new_a_2d[2,0], 2)

        new_a_2d[3,0] = 7
        self.assertEqual(new_a_2d[3,0], 7)

        new_a_2d[3] = np.array([-1, -2, 0])
        self.assertListEqual(list(new_a_2d[3]), [-1, -2, 0])
        self.assertTrue(isinstance(new_a_2d, cctk.OneIndexedArray))

        v = new_a_2d[1]
        self.assertTrue(isinstance(v, np.ndarray))

        new_a_2d[1] = np.asarray([4, 4, 4])
        self.assertListEqual(list(new_a_2d[1]), [4, 4, 4])

    def test_smart_indexing(self):
        a1 = cctk.OneIndexedArray([1,2,3])
        a2 = cctk.OneIndexedArray([4,5,6])

        a3 = np.hstack((a1, a2))
        self.assertListEqual(list(a3), [1, 2, 3, 4, 5, 6])
        self.assertTrue(isinstance(a3, np.ndarray))

        a3 = cctk.OneIndexedArray(a3)
        self.assertListEqual(list(a3), [1, 2, 3, 4, 5, 6])
        self.assertTrue(isinstance(a3, cctk.OneIndexedArray))

        a4 = np.vstack((a3,a3))
        a4 = cctk.OneIndexedArray(a4)
        self.assertTrue(a4.shape[0] == 2)
        self.assertTrue(a4.shape[1] == 6)
        self.assertListEqual(list(a4[1]), [1, 2, 3, 4, 5, 6])

        idx = [True, True, False, True, False, False]
        b1 = a3[idx]
        self.assertListEqual(list(b1), [1, 2, 4])
        self.assertTrue(isinstance(b1, cctk.OneIndexedArray))

        idx = [1, 2, 4]
        b2 = a3[idx]
        self.assertListEqual(list(b1), [1, 2, 4])
        self.assertTrue(isinstance(b1, cctk.OneIndexedArray))

if __name__ == '__main__':
    unittest.main()
