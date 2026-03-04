from enum import Enum
import errno
from importlib.metadata import version
from pathlib import Path
import logging

import numpy as np

logger = logging.getLogger(__name__)

class ExtendedEnum(Enum):

    @classmethod
    def values(cls):
        return [m.value for m in cls]

def ensure_ndarray(v):
    '''
    Return numpy array of v, given array, int/float, list, tuple
    '''
    if isinstance(v, np.ndarray):
        return v

    if isinstance(v, (int, float)):
        return np.array([v])

    return np.array([*v])


def parameter_to_array(parameter, graph_keys):
    '''Parameter argument to numpy array

    Parameters
    ----------
    parameter : array, int, float or dict
        if array:
            make sure length is n
        if int or float:
            return an array of length n with value
        if dict:
            returns an array of length n with values set according to
        keys (nodes indexed by graph_keys), and the rest set to 'default' key

    graph_keys : dict

    Returns
    ----------

    parameter_array: numpy array
        array of length n
    '''

    if isinstance(parameter, np.ndarray) and \
        (len(parameter) == len(graph_keys)):
        return parameter

    parameter_array = np.ones(len(graph_keys))

    if isinstance(parameter, (int, float)):
        parameter_array = parameter_array * parameter

    elif isinstance(parameter, dict):
        if 'default' in parameter.keys():
            parameter_array = parameter_array * parameter['default']
        for key, value in parameter.items():
            if key == 'default':
                continue
            parameter_array[graph_keys[key]] = value
    else:
        logger.warning("'parameters' must be int, float, or dict.")
        parameter_array = parameter_array * 1
    return parameter_array


def file_writable(path):
    '''Checks if path is writeable. If not, attempts to provide reason, and raises
    an Exception.

    Parameters
    ----------
    path : str or Path
        Path to file to check for writability. If the file does not exist, the
        function will check if it can be created.

    Returns
    -------
    None

    '''
    path = Path(path)

    if path.exists():
        logger.warning('%s already exists and will be overwritten!', path)

    try:
        with open(path, 'wb'):
            pass
    except IOError as e:
        if e.errno == errno.EACCES:
            logger.error('No permission to write to %s.', path)
        elif e.errno == errno.EISDIR:
            logger.error('%s is a directory.', path)
        raise e

def get_pkg_version():
    '''Get version of booldog package

    Returns
    -------
    version : str
        Version of booldog package
    '''
    return version('booldog')


