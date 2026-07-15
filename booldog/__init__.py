# Copyright (C) 2020-2026 National Institute of Biology, Slovenia
# Author: Carissa Bleker
# Contact: carissa.bleker@nib.si
"""`booldog`: A Python package for analyses of Boolean and semi-quantitative Boolean networks.

This top-level package module re-exports the public API (`BoolDogModel`,
`BooleanStateSpace`), exposes the installed package version as
`__version__` (via `booldog.utils.misc.get_pkg_version`), and configures
package-wide logging (via `booldog.utils.logger.setup_logger`) on import.
"""

import sys
import logging

from booldog.utils.misc import get_pkg_version
from booldog.utils.logger import setup_logger
from booldog.network import BoolDogModel
from booldog.simulation_result import BooleanStateSpace

assert sys.version_info >= (3, 10)

__version__ = get_pkg_version()

setup_logger()
logger = logging.getLogger(__name__)
logger.info("BoolDog version: %s", __version__)
