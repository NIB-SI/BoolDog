'''Defines `booldog.network.BoolDogModel`, the central class of the package.

`BoolDogModel` carries almost no logic of its own: it is a composition of
mixins (Boolean-network operations, structural modifications, continuous/ODE
simulation, and I/O), each implemented in its own subpackage. This module
only owns `__init__`, node/index bookkeeping, and the `primes` property.
'''

import logging

from pyboolnet import file_exchange

from booldog.classes import BoolDogNode
from booldog.boolean import BooleanNetworkMixin, BooleanNetworkModificationMixin
from booldog.continuous import ContinuousMixin
from booldog.io import BoolDogModelIOFromMixin, BoolDogModelIOToMixin

logger = logging.getLogger(__name__)

# hide runtime warnings (divide by zero, multiply by inf)
# occur in squad ode init
# TODO This is a bad idea...
# warnings.filterwarnings("ignore",
#                         message="divide by zero encountered in true_divide")
# warnings.filterwarnings("ignore",
#                         message="invalid value encountered in multiply")

class BoolDogModel(BooleanNetworkMixin, BooleanNetworkModificationMixin, ContinuousMixin, BoolDogModelIOFromMixin, BoolDogModelIOToMixin):
    '''A class to represent a Boolean network.

    This class is a composition of mixins: `booldog.boolean.BooleanNetworkMixin`
    (read-only Boolean-network operations), `booldog.boolean.BooleanNetworkModificationMixin`
    (structural edits), `booldog.continuous.ContinuousMixin` (ODE-based continuous
    simulation), and `booldog.io.BoolDogModelIOFromMixin`/`booldog.io.BoolDogModelIOToMixin`
    (reading/writing various network exchange formats). This class (`network.py`) itself
    only defines `__init__`, node/index bookkeeping, and the `primes` property.
    '''

    def __init__(self, nodes=None, primes=None, modelinfo=None):
        '''Initialise a Boolean network.

        Parameters
        ----------
        nodes : iterable of BoolDogNode or dict
            Iterable of `booldog.classes.BoolDogNode` objects representing the nodes.
            Each element may also be a dict with (at least) 'identifier' and 'rule'
            keys, which is converted to a `BoolDogNode` automatically. A node with
            no rule (falsy `rule`) is assumed to be an input node, and its rule is
            set to its own identifier.
        primes : dict, optional
            Dictionary of prime implicants, keyed by node identifier. If not given,
            the prime implicants are computed from the node rules (see `get_primes`).
        modelinfo : BoolDogModelInfo, optional
            Model metadata. See `booldog.classes.BoolDogModelInfo` for more
            information. If not given, `self.modelinfo` is not set at all (there
            is no default `BoolDogModelInfo` instance created). **Known bug:**
            some code elsewhere (e.g. `booldog.io.cytoscape`,
            `booldog.simulation_result.continuous_result`) accesses
            `self.modelinfo` unconditionally, so constructing a model
            without `modelinfo` and later hitting one of those code paths
            raises `AttributeError`. See ``KNOWN_BUGS.md``.

        Notes
        -----
        For information on the prime implicants format, see
        `PyBoolNet: prime implicants
        <https://pyboolnet.readthedocs.io/en/latest/Manual.html#prime-implicants>`_.
        '''
        if nodes is None:
            raise ValueError("Nodes must be provided.")

        for node in nodes:
            if not isinstance(node, BoolDogNode):
                # try to make it a BoolDogNode
                try:
                    node = BoolDogNode(**node)
                except Exception as e:
                    raise ValueError(
                        "Nodes must be of type BoolDogNode or dict with keys 'identifier' and 'rule'.") from e

            # if node has no rule, assume it is an 'input' node
            if not node.rule:
                logger.info("Node '%s' has no rule. Assuming 'input' node.", node.identifier)
                node.rule = node.identifier

        self.nodes = {node.identifier: node for node in nodes}
        '''dict of str to BoolDogNode : the network's nodes, keyed by identifier.'''

        if primes is None:
            self._primes = self.get_primes()
        else:
            self._primes = primes
        self._primes_cached = True

        self._set_node_ids_and_index()

        if modelinfo is not None:
            self.modelinfo = modelinfo
            '''BoolDogModelInfo : model metadata, only set if `modelinfo` was provided.'''

        self.modifications = []
        '''list of Modification : audit trail of structural modifications (see
        `booldog.boolean.modifications.Modification`) applied to this network via
        `add_node`/`remove_node(s)`/`update_node`/`modify_network`.'''

        logger.info("Created Network with %i nodes.", self.n)

    def _set_node_ids_and_index(self):
        '''Set the sorted node identifiers and index mapping to integers.

        Called on initialisation, and again after any structural modification
        (see `booldog.boolean.modifications.BooleanNetworkModificationMixin._update_model_object`)
        to keep `node_ids`/`index` consistent with `self.nodes`.
        '''

        # use to maintain consistent node ordering
        self.node_ids = tuple(sorted(self.nodes.keys()))
        '''tuple of str : node identifiers, sorted, defining a consistent node ordering.'''

        # create index mapping (id -> int index)
        self.index = {node_id: i for i, node_id in enumerate(self.node_ids)}
        '''dict of str to int : maps each node identifier to its integer index in
        `node_ids`, used for indexing arrays (e.g. state vectors, ODE parameter arrays).'''

    @property
    def n(self):
        '''Number of nodes in the network'''
        return len(self.nodes)

    @property
    def primes(self):
        '''Prime implicants of the Boolean network. See
        `PyBoolNet:prime implicants
        <https://pyboolnet.readthedocs.io/en/latest/Manual.html#prime-implicants>`_
        for more information.
        '''

        # cache primes to avoid recomputing them every time. This is important for
        # larger networks, where computing primes can take a while.
        if not self._primes_cached:
            self._primes = self.get_primes()
            self._primes_cached = True
        return self._primes

    #############################
    #       REPRESENTATION      #
    #############################

    def __repr__(self):
        '''String representation of the model, showing its class and node count.'''
        return f"{self.__class__} with {self.n} nodes"

    def __len__(self):
        '''Number of nodes in the network (same as ``n``).'''
        return self.n
