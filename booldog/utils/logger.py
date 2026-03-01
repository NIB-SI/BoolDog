import logging
import sys

def setup_logger(level=logging.DEBUG):
    ''' Set up logger for BoolDog
    '''

    # Set up logging
    logger = logging.getLogger('booldog')
    logger.setLevel(level)

    # create console handler and set level to debug
    ch = logging.StreamHandler()
    # ch.setLevel(level)

    # create formatter
    formatter = logging.Formatter('%(levelname)s %(asctime)s %(name)s:%(funcName)s %(message)s')

    # add formatter to ch
    ch.setFormatter(formatter)

    # add ch to logger
    logger.addHandler(ch)

    logger.propagate = False

