from enum import Enum
import errno
from importlib.metadata import version
from pathlib import Path
import logging

import numpy as np

logger = logging.getLogger(__name__)

class ExtendedEnum(Enum):
    '''`enum.Enum` subclass that adds a `values` classmethod for retrieving
    all member values as a plain list.

    Used as the base class for
    `booldog.boolean.modifications.ModificationTypes`.
    '''

    @classmethod
    def values(cls):
        '''Return the values of all members of this enum.

        Returns
        -------
        list
            The ``.value`` of each member of the enum, in definition
            order.
        '''
        return [m.value for m in cls]

def ensure_ndarray(v):
    '''Coerce *v* to a `numpy.ndarray`.

    Parameters
    ----------
    v : numpy.ndarray, int, float, or iterable
        The value to coerce. If *v* is already a `numpy.ndarray` it is
        returned unchanged. If *v* is an `int` or `float`, a length-1
        array containing *v* is returned. Otherwise, *v* is assumed to be
        an iterable (e.g. `list`, `tuple`, generator) and is unpacked into
        a new array via ``np.array([*v])``.

    Returns
    -------
    numpy.ndarray
        *v* coerced to an array, as described above.
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
        Mapping from node identifier to its (0-based) index in the
        returned array, e.g. `booldog.network.BoolDogModel.index`.
        ``len(graph_keys)`` determines *n*, the length of the returned
        array; for a dict *parameter*, its keys are looked up in
        *graph_keys* to find which index to assign each value to.

    Returns
    ----------

    parameter_array: numpy array
        array of length n

    Notes
    -----
    If *parameter* is already a `numpy.ndarray` of length
    ``len(graph_keys)``, it is returned unchanged, bypassing all the cases
    below. Otherwise, an all-ones array of length ``len(graph_keys)`` is
    built and then:

    * if *parameter* is an ``int``/``float``, every entry is set to that
      value;
    * if *parameter* is a ``dict``, every entry is first set to
      ``parameter['default']`` if a ``'default'`` key is present, then
      entries are overwritten per-key from *parameter* (using
      *graph_keys* to map each key to its index); keys of *parameter*
      absent from *graph_keys* raise `KeyError`, and keys of *graph_keys*
      absent from *parameter* keep the all-ones/default value;
    * for any other input -- including a `numpy.ndarray` whose length
      does *not* match ``len(graph_keys)`` -- a warning is logged and an
      all-ones array is returned, silently ignoring the supplied value.
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

    Raises
    ------
    IOError
        If *path* cannot be opened for writing. If the underlying error
        is a permission error (``errno.EACCES``) or *path* is a directory
        (``errno.EISDIR``), a descriptive message is logged first via
        `logging.Logger.error`; the original exception is then re-raised
        unchanged in all cases.

    Notes
    -----
    A warning is logged (but no exception raised) if *path* already
    exists, since it will be overwritten.

    The writability check itself is performed by opening *path* for
    writing (``open(path, 'wb')``) and immediately closing it again
    without writing any data. As a side effect, if *path* did not exist
    it is created as an empty file, and if it did exist its contents are
    truncated. Callers that need the existing content of *path*
    preserved should read/copy it before calling this function.
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

    Notes
    -----
    Reads the version from the installed ``BoolDog`` distribution's
    metadata via `importlib.metadata.version`. This requires the package
    to be installed (e.g. via `pip`/`poetry`) with its distribution
    metadata available; it will raise
    `importlib.metadata.PackageNotFoundError` if ``booldog`` has no
    installed distribution metadata to read.
    '''
    return version('booldog')


