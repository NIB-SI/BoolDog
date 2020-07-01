# https://github.com/hklarner/PyBoolNet/blob/master/Docs/Sphinx/source/Development.rst
import os, sys
#PyBoolNet_path = os.path.join(os.path.abspath("."), "PyBoolNet/")
#sys.path.insert(0, PyBoolNet_path)

from .io import import_graphml
from .io import import_bnet

from .base_class import SquadRegulatoryNetwork