''' Read/write PyBoolNet prime implicants (dict or JSON file).
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
    ''' Parse PyBoolNet prime implicants (from a dict or a JSON file) into
    the data needed to construct a :py:class:`BoolDogModel`.

    The primes are converted to bnet format internally (via
    ``pyboolnet.file_exchange.primes2bnet``) to derive each node's rule
    string; the parsed primes themselves are also kept and returned as-is.

    Parameters
    ----------
    primes_input : str or dict
        Dictionary of primes (PyBoolNet prime-implicant format), or a path
        to an existing file containing primes saved in JSON format (read
        via ``pyboolnet.file_exchange.read_primes``).

    Returns
    -------
    data : dict
        Dictionary with keys ``"nodes"`` (list of :class:`BoolDogNode`,
        one per key in ``primes``), ``"modelinfo"`` (:class:`BoolDogModelInfo`,
        with ``source_format`` set to ``"primes"`` and ``source`` set to
        ``primes_input`` when it is a file path, else ``None``), and
        ``"primes"`` (the parsed prime-implicant dict). Suitable for
        ``BoolDogModel(**data)``.

    Raises
    ------
    ValueError
        If ``primes_input`` is neither a dict nor a path to an existing
        file.
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
    '''Save a model's prime implicants as a formatted JSON file.

    Unlike ``pyboolnet.file_exchange.write_primes``, this always writes
    formatted (``sort_keys=True, indent=2``) JSON and does not support
    returning the JSON as a string; ``outfile`` is required.

    Parameters
    ----------
    network : BoolDogModel
        A BoolDog object representing a Boolean network; its ``primes``
        property is what gets written.
    outfile : str or Path
        File name/path to write primes to. Overwritten if it already
        exists.

    Returns
    -------
    None
    '''
    with open(outfile, "w", encoding="utf-8") as fp:
        json.dump(network.primes, fp, sort_keys=True, indent=2)
