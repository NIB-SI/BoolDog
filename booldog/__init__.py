# Copyright (C) 2020-2021 National Institute of Biology, Slovenia
# Author: Carissa Bleker
# Contact: carissa.bleker@nib.si
"""`booldog`: A Python package for analyses of Boolean and semi-quantitative Boolean networks"""

import sys
from .utils.logger import setup_logger

setup_logger()

assert sys.version_info >= (3, 10)

from .network import BoolDogModel
from booldog.simulation_result import BooleanStateSpace

import logging

logger = logging.getLogger(__name__)

logger.info("BoolDog version: %s", __version__)

__all__ = ['BoolDogModel']

