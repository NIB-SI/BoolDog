'''


See: https://realpython.com/primer-on-python-decorators/
'''

import functools
import time
from booldog.classes import BoolDogNode


# Decorator template:
# def decorator(func):
#     @functools.wraps(func)
#     def wrapper_decorator(*args, **kwargs):
#         # Do something before
#         value = func(*args, **kwargs)
#         # Do something after
#         return value
#     return wrapper_decorator


def timer(func):
    """Print the runtime of the decorated function"""
    @functools.wraps(func)
    def wrapper_timer(*args, **kwargs):
        start_time = time.perf_counter()
        value = func(*args, **kwargs)
        end_time = time.perf_counter()
        run_time = end_time - start_time
        print(f"Finished {func.__name__}() in {run_time:.4f} secs")
        return value
    return wrapper_timer


def debug(func):
    """Print the function signature and return value"""
    @functools.wraps(func)
    def wrapper_debug(*args, **kwargs):
        args_repr = [repr(a) for a in args]
        kwargs_repr = [f"{k}={repr(v)}" for k, v in kwargs.items()]
        signature = ", ".join(args_repr + kwargs_repr)
        print(f"Calling {func.__name__}({signature})")
        value = func(*args, **kwargs)
        print(f"{func.__name__}() returned {repr(value)}")
        return value
    return wrapper_debug

def silence_stdout(func):
    """Silence the standard output of the decorated function."""
    import sys
    import os

    @functools.wraps(func)
    def wrapper_silence_stdout(*args, **kwargs):
        # Save the current stdout so we can restore it later
        original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')
        try:
            value = func(*args, **kwargs)
        finally:
            sys.stdout.close()
            sys.stdout = original_stdout
        return value
    return wrapper_silence_stdout

# decorator that tests if class method argument is a valid node, if a BoolenaNode object,
# changes the argument to be the node.identifier
def validate_node_argument(func):
    '''Decorator to validate node arguments for Boolean network methods.'''

    def wrapper_validate_node_argument(self, node_id, *args, **kwargs):
        if isinstance(node_id, BoolDogNode):
            node_id = node_id.identifier
        if node_id not in self.node_ids:
            raise ValueError(f"{node_id} is not a node identifier of the model.")
        return func(self, node_id, *args, **kwargs)

    return wrapper_validate_node_argument

