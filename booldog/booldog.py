''' TODO doc string'''

import logging, warnings

from booldog.boolean_network import BooleanNetwork
from booldog.semi_quantitative import ContinuousMixin
from booldog.io import WriteMixin
from booldog.io import read
from booldog.ode_factory import ode_factory

logger = logging.getLogger(__name__)

# hide runtime warnings (divide by zero, multiply by inf)
# occur in squad ode init
# TODO This is a bad idea...
# warnings.filterwarnings("ignore",
#                         message="divide by zero encountered in true_divide")
# warnings.filterwarnings("ignore",
#                         message="invalid value encountered in multiply")


class Network(BooleanNetwork, ContinuousMixin, WriteMixin):
    '''A class to represent a Boolean network.

    '''

    def __init__(self, graph=None, **kwargs):
        '''Initialise a Boolean network.

        Parameters
        ----------
        graph : `booldog.BooleanNetwork` or dict or str
            If not a :py:class:`booldog.BooleanNetwork` instance, then correct
            input for initializing a :py:class:`booldog.BooleanNetwork`.

        Other Parameters
        ----------
        **kwargs
            In the case `graph` is not a BooleanNetwork instance, additional
            keyword arguments passed to :py:class:`booldog.BooleanNetwork`.
        '''
        if not isinstance(graph, BooleanNetwork):
            graph = read(graph, **kwargs)

        super().__init__(graph.primes, graph.nodes, graph.index, graph.n)


        self.set_model_source

        logger.info("Created Network with %i nodes", self.n)

    def transform_bool_to_continuous(self,
                                     transform="normalisedhillcube",
                                     **kwargs):
        '''Generate an ODE from RegulatoryNetwork/Boolean graph.

        Note that the BooleanNetwork object is kept in memory as the primes of
        the Boolean network. This means that importing the graph may take a
        while, depending on the size of the network.

        Parameters
        ----------
        transform : str
            One of accepted transforms. See `booldog.ode.transforms` for
            accepted options.

        Other Parameters
        ----------------
        **kwargs
            Additional keyword arguments passed to :py:func:`booldog.ODE`.


        Returns
        ----------


        '''
        ode_system = ode_factory(self, transform, **kwargs)
        return ode_system

    def set_model_source(self, source):
        self.source = source