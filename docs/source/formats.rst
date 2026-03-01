============
Formats
============


.._boolnet_format:
Boolnet
=======

The boolnet (.bnet) format is a simple text format to represent Boolean
networks. Each line contains a target variable and its update function (as a Boolean expression). The symbols ``&``, ``|`` and ``!`` are respectively used for conjunction, disjunction and negation. The target is separated from the update rule by a comma.

The following is an example of a boolnet file:

.. code-block:: bash

    targets, factors
    A, B & (C | D)
    B, !A
    C, A | D & !B

The header of ``targets, factors`` is optional. In addition, comments start with a hash (``#``).

There should only be one line per target variable.

Node names should not contain special characters, including:
"." (period).

Sources
-------
* `BoolNet package vignette <https://rdrr.io/cran/BoolNet/f/inst/doc/BoolNet_package_vignette.Snw.pdf>`_
* `pyboolnet documentation <https://pyboolnet.readthedocs.io/en/master/quickstart.html?highlight=bnet#boolean-networks>`_


Primes
======

A Boolean elementary conjunction `f` is an implicant of a target variable `A` if `f(X) = 1` implies `A(X) = 1`. As an example, the implicants of `A` above includes `B & C`. *Prime* implicants are the shortest of such clauses.

In the vernacular of pyboolnet, 1 implicants correspond to all clauses that imply the expression is true, while 0 implicants correspond to all clauses that are false.

Prime implicants are used as another representation of a Boolean network. They are represented as lists of length two, with the first entry being the 0 prime implicants and the second being the 1 prime implicants. The implicants themselves are each represented by dictionaries, with the key as the component name and the value as 0 or 1, depending whether the component is negated or not.

The previous BooolNet network formatted as primes is as following:

.. code-block:: python

    {
        'A': [
            [ # 0 prime implicants of A
                {'C': 0, 'D': 0},   # 1st 0 prime implicant of A
                {'B': 0}            # 2nd 0 prime implicant of A
            ],
            [ # 1 prime implicants of A
                {'B': 1, 'D': 1},   # 1st 1 prime implicant of A
                {'B': 1, 'C': 1}    # 2nd 1 prime implicant of A
            ]
        ],
        'B': [
            [
                {'A': 1}
            ],
            [
                {'A': 0}
            ]
        ],
        'C': [
            [
                {'A': 0, 'D': 0},
                {'A': 0, 'B': 1}
            ],
            [
                {'B': 0, 'D': 1},
                {'A': 1}
            ]
        ],
        'D': [
            [
                {'D': 0}
            ],
            [
                {'D': 1}
            ]
        ]
    }

Prime implicants can be saved as a JSON file.


Sources
-------
* `pyboolnet documentation <https://pyboolnet.readthedocs.io/en/master/manual.html>`_
* H. Klarner, A. Bockmayr and H. Siebert. (2015). *Computing maximal and minimal trap spaces of Boolean networks.* Natural computing, 14(4). `https://doi.org/10.1007/s11047-015-9520-7 <https://doi.org/10.1007/s11047-015-9520-7>`_
* Crama, Y., & Hammer, P. L. (2011). *Boolean functions: Theory, algorithms, and applications.* Cambridge University Press.



Graphml
=======


Sources
-------
*



SBML-qual
=========


Sources
-------
*



More to come...
===============