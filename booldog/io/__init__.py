''' Mixins that expose ``BoolDogModel.from_XXX``/``to_XXX`` constructors and
exporters, dispatching to the reader/writer functions implemented in the
sibling ``booldog.io`` modules (one module per supported file/object format).
'''
import logging
from io import StringIO

import json

from pyboolnet import file_exchange

from .networkx import booldog2networkx
from .igraph import booldog2igraph
from .cytoscape import booldog2cytoscape

from .bnet import read_bnet, write_bnet

from .interaction_networks import read_sif, read_graphml, read_igraph, read_networkx, read_interactions

from .sbml import read_sbmlqual, write_sbmlqual
from .primes import read_primes, write_primes
from .tabularqual import read_tabularqual

from booldog.utils import file_writable

logger = logging.getLogger(__name__)


class BoolDogModelIOFromMixin:
    """
    Mixin providing BoolDogModel.from_XXX constructors
    for supported Boolean / regulatory network formats.
    """

    # ---- Core helper -------------------------------------------------

    @classmethod
    def _from_reader(cls, reader, *args, **kwargs):
        """
        Generic constructor wrapper around a reader function.

        Reader functions should take ``*args`` and ``**kwargs`` as input and return a
        dictionary of "nodes", "modelinfo" and "primes" (primes and modelinfo
        are optional, but "nodes" is required).
        """

        logger.debug(
            "Creating %s using reader %s",
            cls.__name__,
            reader.__name__,
        )

        data = reader(*args, **kwargs)

        return cls(**data)

    # ---- Boolean network formats --------------------------------------------

    @classmethod
    def from_bnet(cls, *args, **kwargs):
        """Create model from BoolNet .bnet format.

        Forwards all arguments to :func:`booldog.io.bnet.read_bnet`.
        """
        return cls._from_reader(read_bnet, *args, **kwargs)

    @classmethod
    def from_primes(cls, *args, **kwargs):
        """Create model from prime implicants format.

        Forwards all arguments to :func:`booldog.io.primes.read_primes`.
        """
        return cls._from_reader(read_primes, *args, **kwargs)

    @classmethod
    def from_sbmlqual(cls, *args, **kwargs):
        """Create model from SBML-qual file.

        Forwards all arguments to :func:`booldog.io.sbml.read_sbmlqual`.
        """
        return cls._from_reader(read_sbmlqual, *args, **kwargs)

    @classmethod
    def from_tabularqual(cls, *args, **kwargs):
        """Create model from tabular-qual format.

        Forwards all arguments to
        :func:`booldog.io.tabularqual.read_tabularqual`.
        """
        return cls._from_reader(read_tabularqual, *args, **kwargs)

    # ---- Interaction / graph-based formats ----------------------------------

    @classmethod
    def from_interactions(cls, *args, **kwargs):
        """Create model from generic interaction table.

        Forwards all arguments to
        :func:`booldog.io.interaction_networks.read_interactions`.
        """
        return cls._from_reader(read_interactions, *args, **kwargs)

    @classmethod
    def from_sif(cls, *args, **kwargs):
        """Create model from SIF interaction network.

        Forwards all arguments to
        :func:`booldog.io.interaction_networks.read_sif`.
        """
        return cls._from_reader(read_sif, *args, **kwargs)

    @classmethod
    def from_networkx(cls, *args, **kwargs):
        """Create model from networkx graph.

        Forwards all arguments to
        :func:`booldog.io.interaction_networks.read_networkx`.
        """
        return cls._from_reader(read_networkx, *args, **kwargs)

    @classmethod
    def from_igraph(cls, *args, **kwargs):
        """Create model from igraph object or file.

        Forwards all arguments to
        :func:`booldog.io.interaction_networks.read_igraph`.
        """
        return cls._from_reader(read_igraph, *args, **kwargs)

    @classmethod
    def from_graphml(cls, *args, **kwargs):
        """Create model from GraphML file.

        Forwards all arguments to
        :func:`booldog.io.interaction_networks.read_graphml`.
        """
        return cls._from_reader(read_graphml, *args, **kwargs)


class BoolDogModelIOToMixin:
    """
    Mixin providing to_XXX export methods for supported Boolean / regulatory network formats.
    """

    # ---- Core helper --------------------------------------------------------

    @classmethod
    def _to_writer(cls, writer, *args, **kwargs):
        """
        Generic constructor wrapper around a writer function.

        Writer functions should take ``*args`` and ``**kwargs`` as input and return a
        string or write to a file. The first argument should always be the
        model instance (self). An optional keyword argument "outfile" can be
        used to specify an output file path.

        This constructor presents a consistent interface for exporting models
        to various formats, and centralizes any common logic needed.

        Notes
        -----
        None of the ``to_XXX`` methods below currently call this helper; they
        call their writer function directly instead.
        """

        logger.debug(
            "Exporting %s using writer %s",
            cls.__name__,
            writer.__name__,
        )

        if (outfile := kwargs.get("outfile")):
            file_writable(outfile)
            logger.debug("Exporting to file: %s", outfile)
        else:
            logger.debug("Exporting to string output")

        return writer(*args, **kwargs)

    # ---- Boolean network formats --------------------------------------------

    def to_bnet(self, *args, **kwargs):
        """Export model to BoolNet .bnet format."""
        return write_bnet(self, *args, **kwargs)

    def to_primes(self, *args, **kwargs):
        """Export model to prime implicants format."""
        return write_primes(self, *args, **kwargs)

    def to_sbmlqual(self, *args, **kwargs):
        """Export model to SBML-qual file."""
        return write_sbmlqual(self, *args, **kwargs)

    # ---- Object conversions -------------------------------------------------

    def to_networkx(self, *args, **kwargs):
        """Export model to Networkx DiGraph"""
        return booldog2networkx(self, *args, **kwargs)

    def to_igraph(self, *args, **kwargs):
        """Export model to igraph Graph."""
        return booldog2igraph(self, *args, **kwargs)

    # ---- Export to external tools -------------------------------------------

    def to_cytoscape(self, *args, **kwargs):
        """Export model to Cytoscape (via py4cytoscape)."""
        return booldog2cytoscape(self, *args, **kwargs)
