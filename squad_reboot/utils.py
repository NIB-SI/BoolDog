import numpy as np

def ensure_ndarray(v):
    if not type(v) == np.ndarray:
        return np.array([*v])
    else:
        return v


def parameter_to_array(parameter, n, graph_keys, default=1):
    parameter_array = np.ones(n) * default
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
         raise TypeError("'parameters must be int, float, or dict. ")

    return parameter_array


def file_writeable(path):
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














