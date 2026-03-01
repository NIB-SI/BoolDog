'''
Additional functions to parse bnet format

'''
import logging
import os
import re
from io import StringIO

from pyboolnet import file_exchange
from pyboolnet.external.bnet2primes import bnet_file2primes, bnet_text2primes

from booldog.classes import BoolDogNode, BoolDogModelInfo

logger = logging.getLogger(__name__)

##############################
# Utility classes/functions
##############################


BNET_LINE_REGEX = r"^[a-zA-Z]+[a-zA-Z0-9_]*,\s*.+$"
BNET_REGULATORS_REGEX = r"\b[A-Za-z_][A-Za-z0-9_]*\b"

BNET_HEADER = "targets, factors"

class BnetParser:

    def __init__(self, bnet):
        self.bnet = bnet
        self.rules = self._get_rules(bnet)

    def _get_rules(self, bnet_str):
        '''Boolean rules per node.

        '''
        rules = {}

        # nodes that may not have a rule would be missed by just using the
        # targets to collect nodes
        input_nodes = set()

        for line in bnet_str.splitlines():
            line = line.strip()
            if (not line) or (line == BNET_HEADER) or line.startswith("#"):
                continue

            if not re.match(BNET_LINE_REGEX, line):
                logger.error("Bnet line does not conform to bnet format: %s", line)
                logger.debug("Bnet text: %s", bnet_str)
                raise ValueError("Bnet text does not conform to bnet format.")

            logger.debug('Parsing line: %s', line)
            [target, rule] = [s.strip() for s in line.split(",", 1)]
            if target in rules:
                logger.warning("%s already has an update function.", target)
                continue
            rules[target] = rule

            regulators = re.findall(BNET_REGULATORS_REGEX, rule)
            input_nodes.update(regulators)

        for node in input_nodes:
            if node not in rules:
                logger.debug("Adding node %s with no rule.", node)
                rules[node] = ""

        return rules

###############################
# In
###############################

def read_bnet(bnet, node_names=None):
    ''' Generate a BoolDogModel object from a Boolean network in boolnet
    (bnet) format.

    For complete documentation, see :doc:`pyboolnet:modules/file_exchange`.

    Parameters
    ----------
    file : str
        Path to the bnet file

    Returns
    -------
    rn: BoolDog
        An object of type :ref:`py:class:BoolDogModel`.


    Notes
    -----
    The format of the output file is described at :ref:`boolnet_format`.

    '''

    if os.path.exists(bnet):
        with open(bnet, 'r', encoding='utf-8') as f:
            bnet_str = f.read()
        source = bnet
    else:
        bnet_str = bnet
        source = 'object'

    parser = BnetParser(bnet_str)
    nodes = []
    for node_id, rule in parser.rules.items():
        nodes.append(
            BoolDogNode(identifier=node_id,
                        rule=rule,
                        name=node_names.get(node_id, None) if node_names else None))

    modelinfo = BoolDogModelInfo(source=source, source_format="bnet")

    return {
        "nodes": nodes,
        "modelinfo": modelinfo,
        "primes": None
    }

################################
# Out
################################

def write_bnet(model, outfile=None, from_primes=False, header=True, minimize=False):
    ''' Write a BoolDogModel object to a Boolean network in boolnet (bnet) format.

    Parameters
    ----------
    model : BoolDogModel
        A BoolDog object representing a Boolean network.
    outfile : str
        Path to the output file. If None, the output is returned as a string.
    from_primes : bool, default False
        If True, rules are obtained by converting prime implicants to bnet
        format. Otherwise, node rules are written directly.
    header : bool, default True
        If True, include a header line ("target, factors").
    minimize : bool, default False
        If True, minimize rules when converting from primes. Only relevant if
        ``from_primes`` is True.

    Returns
    -------
    str or None
        Returns the bnet string if ``outfile`` is None, otherwise None.

    Notes
    -----
    The output file will be overwritten if it already exists.

    The format of the output file is described at :ref:`boolnet_format`.
    '''
    if from_primes:
        file_exchange.primes2bnet(
            model.primes,
            fname_bnet=outfile,
            header=header,
            minimize=minimize,
        )
        return None

    if outfile is not None:
        f = open(outfile, "w", encoding="utf-8")
        close_file = True
    else:
      f = StringIO()
      close_file = False

    if header:
        f.write(BNET_HEADER + "\n")

    for node_id, node in model.nodes.items():
        if node.rule:
            f.write(f"{node_id}, {node.rule}\n")

    if close_file:
        f.close()
        logger.info("Wrote model as bnet to %s", outfile)
        return None
    else:
        return f.getvalue()
