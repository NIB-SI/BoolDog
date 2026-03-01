import os
import sys
import unittest

import numpy as np

from booldog.utils.misc import ensure_ndarray, parameter_to_array
from booldog.utils.boolean_normal_forms import functions2mindnf

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

    def test_functions2mindnf(self):
        '''Based on pyboolnet tests '''

        functions = {"x": lambda x, y: x or y, "y": lambda x, z: not x and z, "z": lambda x, y, z: sum([x, y, z]) == 1}
        mindnf = functions2mindnf(functions)
        self.assertDictEqual(mindnf, {'x': 'y | x', 'y': 'z&!x', 'z': 'z&!y&!x | !z&!y&x | !z&y&!x'})

    def test_functions2mindnf_depends(self):

        ''' Test the edit '''
        matrix = np.array([[1, 0, 0],
                           [1, 0, 1],
                           [0, 0, 1]])

        n = matrix.shape[0]
        idx = {f"node{x}":x for x in range(n)}

        def function_factory(data, node):
            # actual node dependencies
            args = [f"node{x}" for x in np.nonzero(data[idx[node]])[0]]

            def func(*func_input):
                # derive the network state from the input
                network_state = [0]*n
                for other_node, other_node_state in zip(args, func_input):
                    network_state[idx[other_node]] = other_node_state

                # calculate the node state
                node_state = np.matmul(data[idx[node]], network_state)
                return node_state
            func.depends = args
            return func

        functions = {node:function_factory(matrix,  node) for node in idx.keys()}
        mindnf = functions2mindnf(functions)

        self.assertDictEqual(mindnf, {'node0': 'node0', 'node1': 'node2 | node0', 'node2': 'node2'})

if __name__ == '__main__':
    unittest.main()
