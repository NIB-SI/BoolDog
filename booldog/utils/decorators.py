'''General-purpose function/method decorators used across BoolDog.

See the `Primer on Python Decorators <https://realpython.com/primer-on-python-decorators/>`_
for background on the ``functools.wraps``-based decorator pattern used
throughout this module.
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
    """Print the runtime of the decorated function.

    Parameters
    ----------
    func : callable
        The function to time.

    Returns
    -------
    wrapper_timer : callable
        A wrapped version of *func* that, when called, runs *func*, prints
        its elapsed wall-clock time (via `time.perf_counter`) to stdout,
        and returns *func*'s return value unchanged.
    """
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
    """Print the function signature and return value.

    Parameters
    ----------
    func : callable
        The function to debug.

    Returns
    -------
    wrapper_debug : callable
        A wrapped version of *func* that, when called, prints *func*'s
        name and the ``repr`` of each positional/keyword argument before
        calling it, prints the ``repr`` of its return value afterwards,
        and returns that value unchanged.
    """
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
    """Silence the standard output of the decorated function.

    Redirects `sys.stdout` to `os.devnull` for the duration of the call to
    *func*, restoring the original `sys.stdout` afterwards (even if *func*
    raises).

    Parameters
    ----------
    func : callable
        The function whose stdout output should be suppressed.

    Returns
    -------
    wrapper_silence_stdout : callable
        A wrapped version of *func* that runs it with stdout silenced and
        returns its return value unchanged.
    """
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
    '''Decorator that normalizes and validates the node argument of a
    Boolean-network instance method.

    Intended for methods with signature ``(self, node_id, *args,
    **kwargs)`` on classes that expose a ``self.node_ids`` collection (as
    `booldog.network.BoolDogModel` and its mixins do). Before calling
    *func*, the wrapper:

    * if *node_id* is a `booldog.classes.BoolDogNode`, replaces it with
      its ``identifier`` attribute;
    * checks that the resulting identifier is a member of
      ``self.node_ids``, raising `ValueError` if it is not.

    *func* is then called with the (possibly replaced) plain identifier
    as its second positional argument, plus any remaining
    ``*args``/``**kwargs`` unchanged. This lets decorated methods assume
    *node_id* is always a valid, plain node identifier string, and lets
    callers pass either a `booldog.classes.BoolDogNode` or its identifier
    interchangeably.

    Parameters
    ----------
    func : callable
        The instance method to wrap; must accept ``(self, node_id, *args,
        **kwargs)``.

    Returns
    -------
    wrapper_validate_node_argument : callable
        The wrapped method.

    Raises
    ------
    ValueError
        If the (normalized) *node_id* is not in ``self.node_ids``.
    '''

    def wrapper_validate_node_argument(self, node_id, *args, **kwargs):
        if isinstance(node_id, BoolDogNode):
            node_id = node_id.identifier
        if node_id not in self.node_ids:
            raise ValueError(f"{node_id} is not a node identifier of the model.")
        return func(self, node_id, *args, **kwargs)

    return wrapper_validate_node_argument

