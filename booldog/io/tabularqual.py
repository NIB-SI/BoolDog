'''
Load TabularQual files.

Spec: https://github.com/sys-bio/TabularQual/blob/main/doc/TabularQual_specification_v0.1.2.docx

Uses the ``tabularqual`` package (optional ``sbml`` extra) to convert
TabularQual spreadsheet files to SBML-qual, then loads them using the
SBML-qual reader (:func:`booldog.io.sbml.read_sbmlqual`).
'''

import logging
from tempfile import NamedTemporaryFile

try:
    from tabularqual import convert_spreadsheet_to_sbml
    _TABULARQUAL_AVAILABLE = True
except ImportError as e:
    _TABULARQUAL_AVAILABLE = False

from booldog.io.sbml import read_sbmlqual

logger = logging.getLogger(__name__)

class TabularQualReader:
    '''Reader for TabularQual files, via conversion to SBML-qual.

    Requires the optional ``tabularqual`` package to be installed (part of
    the ``sbml`` extra).

    Parameters
    ----------
    model_path : str
        Path to the TabularQual spreadsheet file.
    '''

    def __init__(self, model_path):
        self.model_path = model_path
        '''str: Path to the TabularQual spreadsheet file, as given at
        construction.'''

    def read(self):
        '''Convert the TabularQual file to SBML-qual (in a temporary file)
        and parse it with :func:`booldog.io.sbml.read_sbmlqual`.

        Returns
        -------
        data : dict
            Dictionary with keys ``"nodes"``, ``"modelinfo"`` and
            ``"primes"`` as returned by
            :func:`booldog.io.sbml.read_sbmlqual`. Suitable for
            ``BoolDogModel(**data)``.

        '''

        with NamedTemporaryFile(delete_on_close=False, suffix=".sbml") as fp:
            convert_spreadsheet_to_sbml(self.model_path, fp.name)
            bn = read_sbmlqual(fp.name)

        return bn

###############################
# In
###############################

def read_tabularqual(model_path):
    '''Parse a TabularQual file into the data needed to construct a
    :py:class:`BoolDogModel`.

    Parameters
    ----------
    model_path : str
        Path to the TabularQual spreadsheet file.

    Returns
    -------
    data : dict
        Dictionary with keys ``"nodes"``, ``"modelinfo"`` and ``"primes"``,
        as returned by :meth:`TabularQualReader.read`. Suitable for
        ``BoolDogModel(**data)``.

    Raises
    ------
    ImportError
        If the optional ``tabularqual`` package is not installed.

    '''

    if not _TABULARQUAL_AVAILABLE:
        raise ImportError("tabularqual is not available.")

    reader = TabularQualReader(model_path)
    return reader.read()
