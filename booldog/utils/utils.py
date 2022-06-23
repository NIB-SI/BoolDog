import numpy as np
from pathlib import Path

def ensure_ndarray(v):
    '''
    Return numpy array of v, given array, int/float, list, tuple
    '''
    if isinstance(v, np.ndarray):
        return v
    elif isinstance(v, (int, float)):
        return np.array([v])
    else:
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
        print("'parameters must be int, float, or dict.")
        parameter_array = parameter_array*1
    return parameter_array


def file_writeable(path):
    '''Checks if path is writable. If not, attempts to print reason, and raises
    an Exception.
    '''
    path = Path(path)

    if path.exists():
        print(f'{path} already exists and will be overwritten')

    try:
        with open(path, 'w') as f:
            pass
    except IOError as x:
        if x.errno == errno.EACCES:
            print(f'No permission to write to {path}')
        elif x.errno == errno.EISDIR:
            print(f'{path} is directory.')
        raise e














