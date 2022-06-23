import os
import sys
import unittest

import numpy as np


sys.path.append("/home/cbleker/research/NIB/squad/BoolDoG")
from booldog.utils.utils import ensure_ndarray, parameter_to_array


class Test(unittest.TestCase):

    def test_ensure_ndarray(self):
        for obj, result in [
            (  1, np.array([1])  ),
            (  [1, 2, 3], np.array([1, 2, 3])  ),
            (  np.array([1, 2, 3]), np.array([1, 2, 3])  )
            ]:
            self.assertTrue(np.array_equal(result, ensure_ndarray(obj), equal_nan=True))

    def test_parameter_to_array(self):
        graph_keys = {
            "A": 0,
            "B": 1,
            "C": 2,
            "D": 3
        }
        for obj, result in [
            (  np.array([1, 2, 3, 4]), np.array([1, 2, 3, 4])  ),
            (  3.3, np.array([3.3, 3.3, 3.3, 3.3])  ),
            (  {"A":1, "B":2, "C":3, "default":5}, np.array([1, 2, 3, 5])  ),
            ]:

            self.assertTrue(np.array_equal(result, parameter_to_array(obj, graph_keys), equal_nan=True))

    def test_file_writeable(self):
        pass

if __name__ == '__main__':
    unittest.main()