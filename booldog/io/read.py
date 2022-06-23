'''Read functions
'''

import os
import sys

import re

import xmltodict
import igraph as ig

import numpy as np
from booldog.utils import *

from collections import defaultdict

from pyboolnet import file_exchange

from pathlib import Path

from booldog.base import RegulatoryNetwork
from booldog.utils import boolean_normal_forms
from booldog.io.squad import SquadInteractions

##############################
#
##############################

def _primes_to_RegulatoryNetwork(primes):

    nodes = tuple(sorted(primes.keys())) # tuple (i.e. immutable)
    index = {i:node for node, i in enumerate(nodes)}
    n = len(nodes)
    return RegulatoryNetwork(primes, nodes, index, n)

def _interactions_to_RegulatoryNetwork(interactions):
    return RegulatoryNetwork(g.primes, g.nodes, g.index, g.n)


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
    rn: :py:class:RegulatoryNetwork
        A RegulatoryNetwork object.

    '''
    if isinstance(primes_input, dict):
        # nothing to do
        primes = primes_input
    else:
        primes = file_exchange.read_primes(primes_input)

    return _primes_to_RegulatoryNetwork(primes)

def read_bnet(bnet):
    ''' Generate a RegulatoryNetwork from a Boolean network in boolnet format.

    For complete documentation, see :doc:`pyboolnet:modules/file_exchange`.

    Parameters
    ----------
    primes_input : str or dict
        Dictionary of primes, of file path to primes saved in JSON format.

    Returns
    -------
    rn: RegulatoryNetwork
        An object of type :ref:`py:class:RegulatoryNetwork`.

    '''

    primes = file_exchange.bnet2primes(bnet)

    return _primes_to_RegulatoryNetwork(primes)


def read_interactions(network):

    interactions = SquadInteractions(network)
    return _interactions_to_RegulatoryNetwork(interactions)


def read_graphml(path,
                 inhibitor_symbol="white_diamond",
                 activator_symbol="standard"):
    '''Reads in a graphml file.


    Parameters
    ----------
    inhibitor_symbol: str, optional
        Symbol of inhibition edges (default="white_diamond")
    activator_symbol: str, optional
        Symbol of activation edges (default="standard")
    '''


    # load graph
    g = ig.Graph.Read_GraphML(path)

    # add edge attributes (i.e. activator or inhibitor)
    with open(path) as f:
        d = xmltodict.parse(f.read())

    D = {}
    N = {}

    for v in d["graphml"]["graph"]["node"]:
        N[v["@id"]] = v["data"]["y:ShapeNode"]["y:NodeLabel"]["#text"]

    for e in d["graphml"]["graph"]["edge"]:
        symbol = e["data"]["y:PolyLineEdge"]['y:Arrows']['@target']
        if symbol == activator_symbol:
            D[e["@id"]] = "+"
        elif symbol == inhibitor_symbol:
            D[e["@id"]] = "-"
        else:
            print(f"Issue with edge ", e["@id"], "arrow symbol \n\t\t{symbol} \n not recognized as activator or inhibitor, perhaps you need to define the 'inhibitor_symbol' ({inhibitor_symbol}) or 'activator_symbol' ({activator_symbol}) in keyword arguments. ")

    for e in g.es():
        e["type"] = D[e["id"]]

    for v in g.vs():
        v["id"] = N[v["id"]]

    pattern = re.compile(r'[^a-zA-Z0-9\-_]')
    search = pattern.search
    for node in g.vs():
        if bool(search(node['id'])):
            print(f"Warning: node names can only contain `{pattern.pattern}`  {node['id']}")
            old_id = node['id']
            node['id'] = re.sub(pattern, '', node['id'])
            #print(f"Warning: renaming node {old_id} --> {node['id']}")

    g_dict = {node["id"]:{} for node in g.vs()}
    for e in g.es():
        source = g.vs()[e.source]["id"]
        target = g.vs()[e.target]["id"]
        sign = e["type"]
        g_dict[source][target] = sign

    interactions = SquadInteractions(g_dict)

    return _interactions_to_RegulatoryNetwork(interactions)