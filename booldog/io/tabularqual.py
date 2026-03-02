'''
Load TabularQual files

Spec: https://github.com/sys-bio/TabularQual/blob/main/doc/TabularQual_specification_v0.1.2.docx

Use tabularqual converter to convert TaabularQual files to SBML-qual, then load them using the SBML-qual reader.
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
    '''Reader for TabularQual files.'''

    def __init__(self, model_path):
        '''Initialise the reader.

        Parameters
        ----------
        model_path : str
            Path to the TabularQual file.

        Returns
        -------
        None

        '''
        self.model_path = model_path

    def read(self):
        '''Read the TabularQual file and return a BoolDogModel.

        Returns
        -------
        BoolDogModel
            The BoolDogModel generated from the TabularQual file.

        '''

        with NamedTemporaryFile(delete_on_close=False, suffix=".sbml") as fp:
            convert_spreadsheet_to_sbml(self.model_path, fp.name)
            bn = read_sbmlqual(fp.name)

        return bn

###############################
# In
###############################

def read_tabularqual(model_path):
    '''Read a TabularQual file and return a BoolDogModel.

    Parameters
    ----------
    model_path : str
        Path to the TabularQual file.

    Returns
    -------
    BoolDogModel
        The BoolDogModel generated from the TabularQual file.

    '''

    if not _TABULARQUAL_AVAILABLE:
        raise ImportError("tabularqual is not available.")

    reader = TabularQualReader(model_path)
    return reader.read()
