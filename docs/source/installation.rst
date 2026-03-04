============
Installation
============

Install from PyPI
=================

To install the latest stable version of BoolDog from PyPI, use the following command:

.. code-block:: bash

   pip install booldog

Optionally, specify a version.


Optional dependencies
=====================

To install optional dependencies, use the following command:

.. code-block:: bash

   pip install booldog[<optional-dependencies>]

Where ``<optional-dependences>`` is a comma-separated list of the following:
* ``sbml``: for SBML support (requires python-libsbml)
* ``networks``: for network analyses (requires networkx and igraph)
* ``graphviz``: for more elaborate STG visualization (requires pygraphviz)

Or:
* ``all``: for all optional dependencies

For details see the `pyproject.toml` file in the repository.

Install from GitHub
===================

Install directly:

.. code-block:: bash

   pip install git+https://github.com/NIB-SI/BoolDog.git

Download and install:

.. code-block:: bash

   git clone https://github.com/NIB-SI/BoolDog.git
   cd BoolDog
   pip install .


Tests
=====

See https://github.com/NIB-SI/BoolDog/tree/main/tests.


Uninstall
=========

.. code-block:: bash

   pip uninstall booldog


Known installation issues
=========================

1. Pyeda not compiling on Windows
---------------------------------
This is an issue with pip and pyeda: https://github.com/cjdrake/pyeda/issues/126.

There is a possible fix, not yet merged: https://github.com/cjdrake/pyeda/pull/153.

To solve, download the appropriate binary from `Christophe Gohlke’s pythonlibs page <https://www.lfd.uci.edu/~gohlke/pythonlibs/#pyeda>`_, and run

.. code-block:: bash

   pip install <downloaded.whl>

Then re-attempt BoolDog installation.
