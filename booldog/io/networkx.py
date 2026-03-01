'''
Function to transform booldog:BoolDogModel to networkx:DiGraph
'''
import logging

import networkx as nx
from pyboolnet.interaction_graphs import primes2igraph
from booldog.io.circuit import booldog2circuit


logger = logging.getLogger(__name__)

def booldog2networkx(model, as_logic_circuit=True):
    '''Export a BoolDog Boolean model to Networkx DiGraph

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
    graph : networkx.DiGraph
        A networkx graph with the same nodes as the input network.
        If `as_logic_circuit` is True, Boolean rules are represented as
        "logical" nodes (and, or, not) and edges.

    Notes
    -----

    See also pyboolnet.interaction_graphs.primes2igraph.

    '''

    if as_logic_circuit:
        g = nx.DiGraph(booldog2circuit(model))
    else:
        g = primes2igraph(model.primes)

        # pyboolnet edges have an attribute "sign" which is a set of 1 (activation)
        # or -1 (inhibition). Based on them add a single "type" attribute
        for _, _, data in g.edges(data=True):
            sign = data.get("sign", set())
            if 1 in sign and -1 in sign:
                data["type"] = "mixed"
            elif 1 in sign:
                data["type"] = "activation"
            elif -1 in sign:
                data["type"] = "inhibition"
            else:
                data["type"] = "none"

        for node_id in g.nodes:
            g.nodes[node_id]["type"] = "species"

    # add node attributes from model (rule, name as label)
    for node_id in g.nodes:
        if node_id not in model.nodes:
            continue
        node = model.nodes.get(node_id, {})
        g.nodes[node_id]["rule"] = node.rule
        g.nodes[node_id]["label"] = node.name

    return g