# Copyright (C) 2020-2021 National Institute of Biology, Slovenia
# Author: Carissa Bleker
# Contact: carissarobyn.bleker@nib.si

"""`booldog`: A Python package for analyses of Boolean and semi-qualitative Boolean networks"""

import sys
import logging


logging.basicConfig(format="%(levelname)s %(asctime)s %(message)s",
                    stream=sys.stdout,
                    level=logging.INFO
)

from .base import RegulatoryNetwork
from .ode import ODE_factory
from . import io

# sphinx stuff
__all__ = ['RegulatoryNetwork', 'ODE_factory']
RegulatoryNetwork.__module__ = "booldog"