'''
This is a modification of pyboolnet.boolean_normal_forms module (originally PyBoolNet.QuineMcCluskey),
specifically the function functions2mindnf, ca line 130.

The original:

    expressions = {}
    for name, func in functions.items():
        inputs = inspect.getargspec(func).args

Modified::

    expressions = {}
    for name, func in functions.items():

        # modification start CB
        if hasattr(func, 'depends'):
            inputs = func.depends
        else:
            inputs = inspect.getargspec(func).args
        # modification end CB


The modification allows the functions in `Functions` argument to accept input
as an k-length list, where k is the number of nodes the node may depend on.

This modification simplifies function creation when using and matrix
multiplication of "Activation" and "Inhibition" matrices when defining node
update functions.

Example:

Suppose we have a matrix that defines the update functions::

    graph_matrix = np.array([[1, 0, 0],
                             [1, 0, 1],
                             [0, 0, 1]])

And the graph state is given by multiplying the matrix  with the current
state::

    node_state = np.matmul(data, network_state)         # line 1

As lambdas (the PyBoolNet format), the functions would be::

    {
    'f1': lambda x1: x1,
    'f2': lambda x1, x2: x1 + x2,
    'f3': lambda x3: x3
    }

I am unaware of a way to automate producing lambdas like these at scale, especially when
`node_state` is composed of more complicated vector operations. The following code
produces functions at scale that have all node states as arguments, calculates
the update, and also identifies the dependencies and make them an attribute of
the function::

    n = data.shape[0]
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

    funcs = {node:function_factory(data,  node) for node in idx.keys()}
'''

# import everything from pyboolnet version
from pyboolnet.boolean_normal_forms import QM

import inspect
import itertools
import logging
from typing import Dict, List

from pyboolnet import NUSMV_KEYWORDS
from pyboolnet.external.bnet2primes import bnet_text2primes

log = logging.getLogger(__name__)

def functions2mindnf(functions: Dict[str, callable]) -> Dict[str, str]:
    """
    Generates and returns a minimal *disjunctive normal form* (DNF) for the Boolean network represented by *functions*.
    The algorithm uses :ref:`Prekas2012 <Prekas2012>`, a Python implementation of the Quine-McCluskey algorithm.

    **arguments**:
        * *functions*: keys are component names and values are Boolean functions

    **returns**:
        * *min_dnf*: keys are component names and values are minimal DNF expressions

    **example**:

        >>> funcs = {"v1": lambda v1,v2: v1 or not v2, "v2": lambda: 1}
        >>> mindnf = functions2primes(funcs)
        >>> mindnf["v1"]
        ((! v2) | v1)
    """

    assert all([inspect.isfunction(f) for f in functions.values()])

    names = functions.keys()

    too_short = [x for x in names if len(x)==1]
    if too_short:
        log.warning(f"variable names must be at least two characters if you want to you NuSMV: names={too_short}")

    forbidden_keywords = [x for x in names if x in NUSMV_KEYWORDS]
    if forbidden_keywords:
        log.warning(f"you are using variable names that are also NuSMV keywords: names={forbidden_keywords}")

    expressions = {}
    for name, func in functions.items():

# >>> BEGIN EDIT
        if hasattr(func, 'depends'):
            inputs = func.depends
        else:
            inputs = inspect.getfullargspec(func).args
# <<< END EDIT

        if not inputs:
            const = func()
            assert const in [0, 1, True, False]
            expressions[name] = "1" if func() else "0"
            continue

        if len(inputs) > 10:
            log.warning(f"computation of prime implicants may take a very long time: name={name}")

        ones, zeros = [], []
        prod = len(inputs) * [[0,1]]
        for i,values in enumerate(itertools.product(*prod)):
            if func(*values):
                ones +=[i]
            else:
                zeros+=[i]

        if not ones:
            expressions[name] = "0"
            continue

        if not zeros:
            expressions[name] = "1"
            continue

        quine = QM(list(reversed(inputs)))
        primes = quine.compute_primes(ones )
        complexity, min_terms = quine.unate_cover(list(primes), ones)
        expressions[name] = quine.get_function(min_terms)

    return expressions

def functions2primes(functions: Dict[str, callable]) -> dict:
    """
    Generates and returns the prime implicants of a Boolean network represented by *functions*.

    **arguments**:
        * *functions*: keys are component names and values are Boolean functions

    **returns**:
        * *primes*: primes implicants

    **example**:

        >>> funcs = {"v1": lambda v1, v2: v1 or not v2,
        ...          "v2": lambda v1, v2: v1 + v2 == 1}
        >>> primes = functions2primes(funcs)
    """

    mindnf = functions2mindnf(functions)
    text = "\n".join([f"{name},\t\t{dnf}" for name, dnf in mindnf.items()])

# >>> BEGIN EDIT
    result = bnet_text2primes(text)
    if result is None:
        raise ValueError("Could not convert functions to primes.")

    return result
# <<< END EDIT











