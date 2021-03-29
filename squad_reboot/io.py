import os
import sys

import xmltodict
import igraph as ig

from PyBoolNet import FileExchange, InteractionGraphs



def import_graphml(path, inhibitor_symbol="white_diamond", activator_symbol="standard"):
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
            print("Issue with edge ", e["@id"], "symbol not activator or inhibitor", symbol)

    for e in g.es():
        e["type"] = D[e["id"]]

    for v in g.vs():
        v["id"] = N[v["id"]]

    g_dict = {node["id"]:{} for node in g.vs()}
    for e in g.es():
        source = g.vs()[e.source]["id"]
        target = g.vs()[e.target]["id"]
        sign = e["type"]
        g_dict[target][source] = sign
    
    return g_dict

def import_bnet(path):
    primes = FileExchange.bnet2primes(path)
    return primes 

def import_primes(path):
    primes = FileExchange.read_primes(path)
    return primes


