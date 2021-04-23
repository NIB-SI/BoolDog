import numpy as np

def ensure_ndarray(v):
    if not type(v) == np.ndarray:
        return np.array([*v])
    else:
        return v


def parameter_to_array(parameter, graph_keys, default=1):
    '''
    Parameter argument to numpy array

    Parameters
    ----------
    parameter : int, float or dict
        if int of float, returns an array of length n with values parameter
        if dict, returns an array of length n with values set according to
        parameter (nodes indexed by graph_keys), and the rest set to default

    graph_keys : dict

	default : int or float, optional
        default value for parameters not set in parameter dict

    Returns
    ----------

    parameter_array: numpy array
        array of length n
    '''

    parameter_array = np.ones(len(graph_keys))
    if isinstance(parameter, (int, float)):
        parameter_array = parameter_array * parameter
    elif isinstance(parameter, dict):
        if default in parameter.keys():
            parameter_array = parameter_array * parameter['default']
        for key, value in parameter.items():
            if key == 'default':
                continue
            parameter_array[graph_keys[key]] = value
    else:
        print("'parameters must be int, float, or dict.")
        parameter_array = parameter_array*default

    return parameter_array


def file_writeable(path):
    '''
    Checks if path is writable. If not, attempts to print reason, and raises
    an Exception.
    '''

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














