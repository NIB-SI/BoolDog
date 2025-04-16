'''Read functions
'''
import logging
from pathlib import Path
from collections import defaultdict

from pyboolnet import file_exchange

from booldog.io.squad import SquadInteractions
from booldog.io.sbml import _SBML_AVAILABLE, sbmlqual2bnet
from booldog.io.graphml import graphml2interactions
from booldog.io.igraph import igraph2interactions
from booldog.io.networkx import networkx2interactions

from booldog.boolean_network import BooleanNetwork

logger = logging.getLogger(__name__)


##############################
#
##############################

forms = [
    "bnet",
    "graphml",
    "interactions",
    "networkx",
    "primes",
    "sbml-qual",
]


def read(
    graph=None,
    format="bnet",
    **kwargs,
):
    '''Generic reader
    '''
    match format:
        case "bnet":
            boolean_graph = read_bnet(graph, **kwargs)
        case "graphml":
            boolean_graph = read_graphml(graph, **kwargs)
        case "sif":
            boolean_graph = read_sif(graph, **kwargs)
        case "interactions":
            boolean_graph = read_interactions(graph, **kwargs)
        case "igraph":
            boolean_graph = read_igraph(graph, **kwargs)
        case "networkx":
            boolean_graph = read_networkx(graph, **kwargs)
        case "primes":
            boolean_graph = read_primes(graph, **kwargs)
        case "sbml-qual":
            boolean_graph = read_sbmlqual(graph, **kwargs)
        case _:
            raise ValueError(f"`format` argument must be of {', '.join(forms)}")

    return boolean_graph

def _primes_to_BooleanNetwork(primes):
    nodes = tuple(sorted(primes.keys()))  # tuple (i.e. immutable)
    index = {i: node for node, i in enumerate(nodes)}
    n = len(nodes)
    return BooleanNetwork(primes=primes, nodes=nodes, index=index, n=n)

def _interactions_to_BooleanNetwork(g):
    return BooleanNetwork(primes=g.primes, nodes=g.nodes, index=g.index, n=g.n)


##############################
#      READ  FUNCTIONS       #
##############################


def read_primes(primes_input):
    ''' Read primes from a dictionary or a json file

    Parameters
    ----------
    primes_input : str or dict
        Dictionary of primes, or a file path to primes saved in JSON format.

    Returns
    -------
    rn: :py:class:BooleanNetwork
        A BoolDog BooleanNetwork object.

    '''
    if isinstance(primes_input, dict):
        # nothing to do
        primes = primes_input
    else:
        primes = file_exchange.read_primes(primes_input)

    return _primes_to_BooleanNetwork(primes)


def read_bnet(file):
    ''' Generate a BoolDog BooleanNetwork object from a Boolean network in boolnet format.

    For complete documentation, see :doc:`pyboolnet:modules/file_exchange`.

    Parameters
    ----------
    primes_input : str or dict
        Dictionary of primes, of file path to primes saved in JSON format.

    Returns
    -------
    rn: BoolDog
        An object of type :ref:`py:class:BoolDog`.

    '''

    primes = file_exchange.bnet2primes(file)

    return _primes_to_BooleanNetwork(primes)


def read_sif(
    file,
    delim="\t",
    header=True,
    **kwargs
):
    '''Reads in a SIF file of interactions

    Parameters
    ----------
    inhibitor_symbol: str, optional
        Symbol of inhibition edges (default="-1")
    activator_symbol: str, optional
        Symbol of activation edges (default="1")



    Note
    ----
    Uses SQUAD logic to obtain Boolean graph
    '''

    interactions = defaultdict(dict)
    with open(file, "r") as handle:
        if header:
            handle.readline()
        for line in handle:
            source, target, interaction = line.rstrip().split(delim)[:3]
            interactions[source][target] = interaction

    interactions_graph = SquadInteractions(interactions, **kwargs)
    return _interactions_to_BooleanNetwork(interactions_graph)

def read_interactions(interactions, **kwargs):
    ''' Create BooleanNetwork from a dictionary of interactions
    '''
    interactions_graph = SquadInteractions(interactions, **kwargs)
    return _interactions_to_BooleanNetwork(interactions_graph)

def read_sbmlqual(file):
    ''' Create BooleanNetwork from a SBML-qual file
    '''
    if _SBML_AVAILABLE:
        bnet = sbmlqual2bnet(file)
        return read_bnet(bnet)

    raise ImportError(
        'libsbml (https://sbml.org/software/libsbml/libsbml-docs/api/python/) '
        'is needed to read models in SBML format. '
        'We suggest you install it using pip. '
    )

def read_graphml(file, **kwargs):
    interactions = graphml2interactions(file, **kwargs)
    interactions_graph = SquadInteractions(interactions, **kwargs)
    return _interactions_to_BooleanNetwork(interactions_graph)

def read_igraph(g, **kwargs):
    ''' Create BooleanNetwork from a igraph object
    '''
    interactions = igraph2interactions(g, **kwargs)
    interactions_graph = SquadInteractions(interactions, **kwargs)
    return _interactions_to_BooleanNetwork(interactions_graph)

def read_networkx(g, **kwargs):
    ''' Create BooleanNetwork from a igraph object
    '''
    interactions = networkx2interactions(g, **kwargs)
    interactions_graph = SquadInteractions(interactions, **kwargs)
    return _interactions_to_BooleanNetwork(interactions_graph)
