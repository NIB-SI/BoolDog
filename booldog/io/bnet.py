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


def _is_bnet_header(line):
    '''Check whether a line is a bnet header line.

    Mirrors BoolNet's own tolerant ``loadNetwork()`` check (per-field,
    whitespace- and case-insensitive), rather than an exact match against
    :data:`BNET_HEADER`: real bnet files vary in header spacing (e.g.
    ``"targets,factors"`` with no space, as used throughout the
    biodivine-boolean-models repository, versus BoolNet's own documented
    ``"targets, factors"``). Only the first two fields are checked;
    BoolNet also allows an optional third ``probabilities`` field for
    probabilistic Boolean networks, which is simply ignored here since
    BoolDog does not otherwise support parsing rules for those.

    Parameters
    ----------
    line : str
        A single, already-stripped line of bnet text.

    Returns
    -------
    bool
    '''
    fields = [f.strip().lower() for f in line.split(",")]
    return len(fields) >= 2 and fields[0] == "targets" and fields[1] in ("functions", "factors")


class BnetParser:
    '''Parse the body of a bnet-format string into per-node Boolean rules.

    Parameters
    ----------
    bnet : str
        Text in bnet format (one ``target, rule`` line per node; blank
        lines, a header line such as ``"targets, factors"`` (see
        :func:`_is_bnet_header` for the accepted variants), and lines
        starting with ``#`` are ignored).
    '''

    def __init__(self, bnet):
        self.bnet = bnet
        '''str: The raw bnet text passed in at construction.'''
        self.rules = self._get_rules(bnet)
        '''dict: Mapping of node identifier to its bnet-format rule string
        (see :meth:`_get_rules`).'''

    def _get_rules(self, bnet_str):
        '''Parse bnet text into a dict of node identifier to rule string.

        Any regulator referenced in a rule that does not itself appear as a
        target elsewhere in the text is also added to the returned dict,
        with an empty string as its rule (i.e. it is treated as an input
        node with no update function).

        Parameters
        ----------
        bnet_str : str
            Text in bnet format.

        Returns
        -------
        rules : dict
            Mapping of node identifier (str) to rule (str). Nodes with no
            defined rule (regulator-only/input nodes) map to ``""``.

        Raises
        ------
        ValueError
            If a non-blank, non-header, non-comment line does not match the
            expected ``target, rule`` bnet line format.
        '''
        rules = {}

        # nodes that may not have a rule would be missed by just using the
        # targets to collect nodes
        input_nodes = set()

        for line in bnet_str.splitlines():
            line = line.strip()
            if (not line) or _is_bnet_header(line) or line.startswith("#"):
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
    ''' Parse a Boolean network in BoolNet (bnet) format into the data
    needed to construct a :py:class:`BoolDogModel`.

    For complete documentation of the bnet format, see
    :doc:`pyboolnet:modules/file_exchange`.

    Parameters
    ----------
    bnet : str
        Either a path to a bnet file, or a string already containing the
        bnet-format text (checked with ``os.path.exists``; if it isn't an
        existing path, it is parsed directly as bnet text).
    node_names : dict, optional
        Mapping of node identifier to a display name, used to populate
        :attr:`BoolDogNode.name` for each parsed node. Nodes not present in
        this dict (or if ``node_names`` is None) get no explicit name (see
        :class:`BoolDogNode`, whose ``name`` defaults to the identifier).

    Returns
    -------
    data : dict
        Dictionary with keys ``"nodes"`` (list of :class:`BoolDogNode`),
        ``"modelinfo"`` (:class:`BoolDogModelInfo`, with ``source_format``
        set to ``"bnet"``), and ``"primes"`` (``None``, since primes are
        not computed by this reader). Suitable for ``BoolDogModel(**data)``.

    Notes
    -----
    The format of the output file is described at
    :doc:`pyboolnet:modules/file_exchange`.

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
        Returns the bnet string if ``outfile`` is None, otherwise writes it
        to ``outfile`` and returns None. This holds for both
        ``from_primes`` values.

    Notes
    -----
    The output file will be overwritten if it already exists.

    The format of the output file is described at
    :doc:`pyboolnet:modules/file_exchange`.
    '''
    if from_primes:
        text_bnet = file_exchange.primes2bnet(
            model.primes,
            fname_bnet=outfile,
            header=header,
            minimize=minimize,
        )
        return text_bnet if outfile is None else None

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
