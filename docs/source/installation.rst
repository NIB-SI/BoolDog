============
Installation
============

Install
=======

Install directly:

.. code-block:: bash

   pip install git+https://github.com/NIB-SI/BoolDog.git

Download and install:

.. code-block:: bash

   git clone https://github.com/NIB-SI/BoolDog.git
   cd BoolDog
   pip install .


Tests
-----

See `https://github.com/NIB-SI/BoolDog/tree/master/tests`.


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

Then re-attempted BoolDog installation.

Dependencies
============

Required dependencies will be automatically installed when installing BoolDog.
The following dependencies are required:

* numpy
* xmltodict
* scipy
* python-igraph
* matplotlib
* pygraphviz (optional)
* pyboolnet
* networkx

