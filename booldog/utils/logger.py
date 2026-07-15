import logging
import sys

def setup_logger(level=logging.INFO):
    '''Configure the ``'booldog'`` logger with a console handler.

    Retrieves the logger named ``'booldog'`` (the parent logger for all
    `booldog` submodules, which each call ``logging.getLogger(__name__)``),
    sets its level to *level*, attaches a `logging.StreamHandler` writing
    to the console with the format
    ``'%(levelname)s %(asctime)s %(name)s:%(funcName)s %(message)s'``, and
    disables propagation to the root logger (``logger.propagate = False``)
    so messages are not also emitted by any handlers attached higher up.

    Parameters
    ----------
    level : int, default logging.INFO
        The logging level to set on the ``'booldog'`` logger (e.g.
        `logging.DEBUG`, `logging.INFO`, `logging.WARNING`).

    Returns
    -------
    None

    Notes
    -----
    Calling this function more than once attaches an additional
    `logging.StreamHandler` to the ``'booldog'`` logger each time, since
    it does not check for or remove handlers that were added by a
    previous call; this results in duplicated console output if
    `setup_logger` is called repeatedly.
    '''

    # Set up logging
    logger = logging.getLogger('booldog')
    logger.setLevel(level)

    # create console handler
    ch = logging.StreamHandler()

    # create formatter
    formatter = logging.Formatter('%(levelname)s %(asctime)s %(name)s:%(funcName)s %(message)s')

    # add formatter to ch
    ch.setFormatter(formatter)

    # add ch to logger
    logger.addHandler(ch)

    logger.propagate = False

    return
