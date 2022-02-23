============
File formats
============

Boolnet
=======


The boolnet (.bnet) format is a simple text format to represent Boolean
networks. Each line contains a target variable and its update function. ``&``, ``|`` and ``!`` are respectively used for conjunction, disjunction and negation. The target is separated from the update rule by a comma.

The following is an example of a boolnet file:

.. code-block:: bash

    A, B & (C | D)
    B, !A
    C, A | D & !B

An optional header of ``targets, factors`` is allowed. In addition, comments start with a hash (``#``).

There should only be one line per target variable.

Sources
-------

* `BoolNet package vignette <https://rdrr.io/cran/BoolNet/f/inst/doc/BoolNet_package_vignette.Snw.pdf>`_
* `pyboolnet documentation <https://pyboolnet.readthedocs.io/en/master/quickstart.html?highlight=bnet#boolean-networks>`_


More to come...
===============