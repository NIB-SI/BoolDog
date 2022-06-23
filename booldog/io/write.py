import os
import sys

import re

import xmltodict
import igraph as ig

import numpy as np

from collections import defaultdict

from pyboolnet import file_exchange

from pathlib import Path

from booldog.utils import *
from booldog.utils import boolean_normal_forms


##############################
#           WRITE            #
##############################


class WriteMixin():

    def write_bnet(self, outfile, minimize=False):
        file_writable(outfile)
        file_exchange.primes2bnet(self.primes, outfile, minimize=True)


    def write_primes(self, outfile):
        file_writable(outfile)
        file_exchange.write_primes(self.primes, outfile)


    def write_sif():


