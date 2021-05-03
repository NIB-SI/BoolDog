# Copyright (C) 2020-2021 National Institute of Biology, Slovenia
# Author: Carissa Bleker
# Contact: carissarobyn.bleker@nib.si

"""`booldog`: A Python package for analyses of Boolean and semi-qualitative Boolean networks"""

from .base import RegulatoryNetwork
from .ode import ODE_factory

# sphinx stuff
__all__ = ['RegulatoryNetwork', 'ODE_factory']
RegulatoryNetwork.__module__ = "booldog"