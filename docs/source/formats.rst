============
Formats
============

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

The header of ``targets, factors`` is optional. If present, it is recognised
regardless of whitespace or case (e.g. ``targets,factors`` also works), with
``functions`` accepted in place of ``factors``, and an optional third
``probabilities`` field (used by BoolNet for probabilistic Boolean networks,
which BoolDog does not otherwise support), matching BoolNet's own tolerant
handling of this header. In addition, comments start with a hash (``#``).

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

SBML-qual and TabularQual
=========================

The SBML-qual format is a standard format (in XML) for representing qualitative models,
including Boolean networks. It is an extension of the Systems Biology Markup
Language (SBML) and allows for the representation of entities, interactions,
and logical rules governing the behaviour of the system.

The TabularQual format is a more user-friendly way to represent the same information as SBML-qual, using a spreadsheet format.

TabularQual is supported by internally converting the TabularQaul spreadsheet to an SBML-qual file,
and parsing the SBML-qual file to extract the relevant information.

Sources
-------
* Chaouiya, C., Bérenguier, D., Keating, S. M., Naldi, A., Van Iersel, M. P., Rodriguez, N., ... & Helikar, T. (2013). SBML qualitative models: a model representation format and infrastructure to foster interactions between qualitative modelling formalisms and tools. BMC systems biology, 7(1), 135. `https://doi.org/10.1186/1752-0509-7-135 <https://doi.org/10.1186/1752-0509-7-135>`_
* `SBML specification <https://sbml.org/documents/specifications/>`_
* `SBML-qual specification <https://sbml.org/documents/specifications/level-3/version-1/qual/>`_

Interactions
============

Pairwise interactions between entities can be imported from:

* Graphml files
* SIF files
* NetworkX DiGraph objects
* igraph Graph objects
* List of interactions (tuples of source, target, and sign)

These can be used to construct a Boolean network by assigning update rules to
each target variable. By default the rules follow the SQUAD convention
(:class:`booldog.io.interaction_logic.SquadLogic`), where the rule of a node
depends on which kinds of regulators it has:

* **no regulators**: constant ``0``
* **activators only**: OR of the activators, e.g. ``A | B``
* **inhibitors only**: AND of the negated inhibitors, e.g. ``!C & !D``
  (the node is on unless an inhibitor is active)
* **both**: the activators OR-ed, AND the negated inhibitors AND-ed, e.g.
  ``(A | B) & (!C & !D)``

This logic can be replaced by passing a custom
:class:`booldog.io.interaction_logic.LogicBuilder` via the ``logic`` argument,
which is accepted by all of the import functions above
(``from_graphml``, ``from_sif``, ``from_networkx``, ``from_igraph`` and
``from_interactions``).

Signs
-----

Every interaction must be signed, i.e. be either an activation or an
inhibition (set via ``activator_symbol`` / ``inhibitor_symbol``); edges of
unknown or dual monotonicity are not supported. Interactions whose sign is
missing or not recognised are handled as follows:

* **Unrecognised sign** (a value matching neither ``activator_symbol`` nor
  ``inhibitor_symbol``, including missing values that igraph fills in as
  ``None``/``NaN``): the interaction is dropped and a warning is logged.
  This applies to all entry points.
* **Missing sign attribute**: an error is raised when the sign cannot be read
  at all, namely when the ``edge_type_key`` attribute is absent from a
  NetworkX edge (``KeyError``), absent from all edges of an igraph/GraphML
  graph (``KeyError``), when a yEd GraphML edge has no arrow head
  (``ValueError``, with ``yEd_arrow_head=True``), or when a SIF line has too
  few columns (``IndexError``).

Note that the sign values are compared as-is: SIF values are always strings
(defaults ``"1"``/``"-1"``), and GraphML values are typed by the ``attr.type``
of their key (e.g. a ``string``-typed key requires
``activator_symbol="1"``, ``inhibitor_symbol="-1"``).

Duplicate interactions
----------------------

If the same (source, target) pair occurs more than once, only the last
occurrence (with a recognised sign) is kept and a warning is logged, stating
whether the signs conflict. In particular, an edge that is both activating
and inhibiting cannot be represented.
