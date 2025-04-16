import os
import sys
import json
import re

from collections import defaultdict

from pyboolnet import file_exchange
from pyboolnet.interaction_graphs import primes2igraph

from pathlib import Path

from booldog.utils import file_writable
from booldog.utils import boolean_normal_forms
from .circuit import booldog2networkx
from .sbml import _SBML_AVAILABLE, booldog2sbmlqual
from .cytoscape import _CYTOSCAPE_AVAILABLE, booldog2cytoscape

##############################
#           WRITE            #
##############################


class WriteMixin():
    '''Mixin class for serialising network object '''

    def write_bnet(self, outfile, **kwargs):

        # IDEA: Keep the original (and updated) logic rules as
        # attribute, and have an option to use them to create the rules
        # transitions, instead of using the bnet/primes DNF

        file_writable(outfile)
        file_exchange.primes2bnet(self.primes, outfile, **kwargs)

    def write_primes(self, outfile):
        file_writable(outfile)
        file_exchange.write_primes(self.primes, outfile)

    def write_sif(self, hypergraph=False):
        g = self.to_networkx(hypergraph=hypergraph)

        # TODO

    def write_sbml(self, outfile, **kwargs):
        ''' '''

        # IDEA: Keep the original (and updated) logic rules as
        # attribute, and have an option to use them to create the sbml -qual
        # transitions, instead of using the bnet/primes DNF

        if _SBML_AVAILABLE:
            return booldog2sbmlqual(self, outfile, **kwargs)

        raise ImportError(
            'libsbml (https://sbml.org/software/libsbml/libsbml-docs/api/python/) '
            'is needed to write models to SBML format. '
            'We suggest you install it using pip. '
        )

    def write_primes2json(self, outfile):
        '''Save primes as formatted JSON file.
        See also pyboolnet.file_exchange.write_primes.

        Parameters
        ==========
        outfile : Path
            File name/path to write primes to.

        '''
        with open(outfile, "w", encoding="utf-8") as fp:
            json.dump(self.primes, fp, sort_keys=True, indent=2)

    def to_networkx(self, as_logic_circuit=False):
        '''Export graph to Networkx DiGraph

        See also pyboolnet.interaction_graphs.primes2igraph.


        '''

        if as_logic_circuit:
            return booldog2networkx(self)

        else:
            return primes2igraph(self.primes)

    def to_igraph(self, hypergraph=False):
        import igraph as ig
        #TODO
        pass

    def to_bnet(self, **kwargs):
        return file_exchange.primes2bnet(self.primes, **kwargs)

    def to_cytoscape(self, **kwargs):
        if  _CYTOSCAPE_AVAILABLE:
            return booldog2cytoscape(self, **kwargs)

        raise ImportError(
            'py4cytoscape (https://py4cytoscape.readthedocs.io/) '
            'is needed to interact with Cytoscape. '
            'We suggest you install it using pip. '
        )
