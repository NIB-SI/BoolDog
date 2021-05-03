============
Installation
============

Install
=======

Download and install:

.. code-block:: bash

   git clone https://github.com/NIB-SI/BoolDoG.git
   cd BoolDoG
   pip install .


Install directly:

.. code-block:: bash

   pip install git+https://github.com/NIB-SI/sBoolDoG.git


Uninstall
=========

.. code-block:: bash

   pip uninstall booldog


Known issues
============

1. Pyeda not compiling on Windows
---------------------------------
This is an issue with pip and pyeda: https://github.com/cjdrake/pyeda/issues/126.

There is a possible fix, not yet merged: https://github.com/cjdrake/pyeda/pull/153. 

To solve, download the appropriate binary from `Christophe Gohlke’s pythonlibs page <https://www.lfd.uci.edu/~gohlke/pythonlibs/#pyeda>`_, and run 

.. code-block:: bash

   pip install <downladed.whl>

Then re-attempted BoolDoG installation. 

Dependencies
============

The following dependencies are required:

* numpy
* xmltodict
* scipy
* python-igraph
* matplotlib
* pygraphviz (optional)
* PyBoolNet
* networkx
