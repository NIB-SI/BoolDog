import numpy as np

def ensure_ndarray(v):
    if not type(v) == np.ndarray:
        return np.array([*v])
    else:
        return v