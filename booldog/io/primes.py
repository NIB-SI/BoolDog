'''
'''
import json
import os

from pyboolnet import file_exchange

from booldog.utils import file_writable
from booldog.classes import BoolDogNode, BoolDogModelInfo

from booldog.io.bnet import BnetParser

##############################
# Utility classes/functions
##############################


###############################
# In
###############################

def read_primes(primes_input):
    ''' Read primes from a dictionary or a json file

    Parameters
    ----------
    primes_input : str or dict
        Dictionary of primes, or a file path to primes saved in JSON format.

    Returns
    -------
    rn: :py:class:BoolDogModel
        A BoolDog BoolDogModel object.

    '''
    if isinstance(primes_input, dict):
        primes = primes_input
        source = None
    elif isinstance(primes_input, str) and os.path.isfile(primes_input):
        primes = file_exchange.read_primes(primes_input)
        source = primes_input
    else:
        raise ValueError(
            "Primes must be a dictionary or a path to a JSON file.")

    bnet = file_exchange.primes2bnet(primes)

    parser = BnetParser(bnet)
    nodes = []
    for node_id, rule in parser.rules.items():
        nodes.append(BoolDogNode(identifier=node_id, rule=rule))

    modelinfo = BoolDogModelInfo(source=source, source_format="primes")

    return {
        "nodes": nodes,
        "modelinfo": modelinfo,
        "primes": primes
    }

################################
# Out
################################

def write_primes(network, outfile):
    '''Save primes as formatted JSON file.
    See also pyboolnet.file_exchange.write_primes.

    Parameters
    ==========
    outfile : Path
        File name/path to write primes to.

    '''
    with open(outfile, "w", encoding="utf-8") as fp:
        json.dump(network.primes, fp, sort_keys=True, indent=2)
