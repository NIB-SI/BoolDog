======================
Miscellaneous settings
======================

Adjust logging level
====================

The default logging level for the application is set to INFO. If you want to see less or more detailed logs, you can change the logging level using the `logging` module from the Python standard library, e.g.:

.. code-block:: python

    import logging
    from booldog import logger
    logger.setLevel(logging.DEBUG)

Sources:
--------
* `logging module documentation <https://docs.python.org/3/library/logging.html>`_



