# Copyright (C) 2020-2021 National Institute of Biology, Slovenia
# Author: Carissa Bleker
# Contact: carissa.bleker@nib.si
"""`booldog`: A Python package for analyses of Boolean and semi-qualitative Boolean networks"""

import logging
import sys

assert sys.version_info >= (3, 10)

logging.basicConfig(
    format="%(levelname)s %(asctime)s %(name)s:%(funcName)s %(message)s",
                    stream=sys.stdout,
                    level=logging.DEBUG)

from . import io
from .booldog import Network

silent_loggers = []
for key in logging.Logger.manager.loggerDict:
    if not key.startswith('booldog'):
        silent_loggers.append(key)
for key in silent_loggers:
    logging.getLogger(key).setLevel(logging.CRITICAL)



# sphinx stuff
__all__ = ['Network', 'ODE_factory']
Network.__module__ = "booldog"

