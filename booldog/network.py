''' TODO doc string'''

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

    Attributes
    ----------
    n : int
        The number of nodes/variables in the network

    primes : dict
        Prime implicants of the Boolean network. See
        `PyBoolNet:prime implicants
        <https://pyboolnet.readthedocs.io/en/latest/Manual.html#prime-implicants>`_
        for more information.

    nodes : tuple of str
        List of `booldog.network.BoolDogNode` objects representing the nodes

    index : dict
        Dictionary of node name to integer index for indexing arrays

    '''

    def __init__(self, nodes=None, primes=None, modelinfo=None):
        '''Initialise a Boolean network.

        Parameters
        ----------
        nodes : iterable of BoolDogNode
            Iterable of `booldog.classes.BoolDogNode` objects representing the nodes.
        primes : dict
            Dictionary of prime implicants. The keys are the node identifiers.
        modelinfo : dict
            Dictionary of model metadata. See `py:class:BoolDogModelInfo` for more information.

        Returns
        -------


        Notes
        -----

        For information on the prime implicants format, see
        `pyboolnet:prime implicants  https://pyboolnet.readthedocs.io/en/master/manual.html#prime-implicants`.
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

        if primes is None:
            self._primes = self.get_primes()
        else:
            self._primes = primes
        self._primes_cached = True

        self._set_node_ids_and_index()

        if modelinfo is not None:
            self.modelinfo = modelinfo

        self.modifications = []

        logger.info("Created Network with %i nodes.", self.n)

    def _set_node_ids_and_index(self):
        '''Set the sorted node identifiers and index mapping to integers.'''

        # use to maintain consistent node ordering
        self.node_ids = tuple(sorted(self.nodes.keys()))

        # create index mapping (id -> int index)
        self.index = {node_id: i for i, node_id in enumerate(self.node_ids)}

    @property
    def n(self):
        '''Number of nodes in the network'''
        return len(self.nodes)

    @property
    def primes(self):
        '''Prime implicants of the Boolean network. See
        `PyBoolNet:prime implicants` <https://pyboolnet.readthedocs.io/en/latest/Manual.html
        #prime-implicants> for more information.
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
        return f"{self.__class__} with {self.n} nodes"

    def __len__(self):
        return self.n
