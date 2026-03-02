'''
Function to transform booldog:Network to DiGraph
'''
import logging
from pyboolnet.interaction_graphs import primes2igraph

try:
    import igraph as ig
    _IGRAPH_AVAILABLE = True
except ImportError as e:
    _IGRAPH_AVAILABLE = False

from booldog.io.circuit import booldog2circuit

logger = logging.getLogger(__name__)

def booldog2igraph(model, as_logic_circuit=True):
    '''Export a BoolDog Boolean model to an igraph graph object.

    Parameters
    ----------
    model : booldog:BoolDogModel
        A BoolDog object representing a Boolean network.

    as_logic_circuit: bool
        If True, the graph is exported as a logic circuit (Boolean rules
        are represented as "logical" nodes (and, or, not) and edges.
        Otherwise, it is exported as a directed interaction graph. Default is False.

    Returns
    -------
    graph : igraph.Graph
        A igraph.Graph object with the same nodes as the input network.
        If `as_logic_circuit` is True, Boolean rules are represented as
        "logical" nodes (and, or, not) and edges.

    Notes
    -----
    See also pyboolnet.interaction_graphs.primes2igraph.

    Implemented via conversion to Networkx DiGraph, then to igraph Graph,
    to reuse pyboolet function and logic circuits code.
    '''

    if not _IGRAPH_AVAILABLE:
        raise ImportError("igraph is not available.")

    if as_logic_circuit:
        g = booldog2circuit(model)
    else:
        g = model.to_networkx(as_logic_circuit=False)

    return ig.Graph.from_networkx(g)
