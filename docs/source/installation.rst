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
* ``biomodels``: for fetching models from the BioModels database (requires requests)

Or:

* ``all``: for all optional dependencies

For example, to install all optional dependencies:

.. code-block:: bash

   pip install booldog[all]

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

1. pyboolnet incompatible with Apple Silicon (ARM) Macs
--------------------------------------------------------
BoolDog depends on `pyboolnet <https://github.com/hklarner/pyboolnet>`_ for
Boolean network analysis. pyboolnet bundles precompiled, Intel-only
binaries (``BNetToPrime``, ``clasp``, ``gringo``, ...), which require
`Rosetta 2 <https://support.apple.com/en-us/102527>`_ to run on Apple
Silicon (M-series) Macs. This has been confirmed via CI on GitHub's Apple
Silicon macOS runners
(`example run <https://github.com/carissableker/pyboolnet/actions/runs/28822384806>`_),
where Rosetta 2 is present and ``BNetToPrime`` instantiation works
correctly.

If you see ``OSError: [Errno 86] Bad CPU type in executable`` as soon as
you instantiate a ``BoolDogModel``, Rosetta 2 is most likely not installed
on your Mac. Per Apple's documentation linked above, it isn't
preinstalled, and normally only installs itself the first time you
*interactively* open an Intel-based app and accept the install prompt,
which never fires for a binary invoked via Python's ``subprocess``.
Installing it manually from a Terminal should resolve this:

.. code-block:: bash

   softwareupdate --install-rosetta --agree-to-license

Separately, even with Rosetta 2 installed, piping two of these translated
binaries together (``gringo | clasp``, used for trap-space computation,
i.e. ``BoolDogModel.steady_states()``) can fail with a "bad input stream"
parse error. In CI, this reproduces reliably on the Apple Silicon macOS
runner and not on Linux or Windows runners, which run the same binaries
natively
(`CI run <https://github.com/carissableker/pyboolnet/actions/runs/28822384806>`_).
The exact mechanism is unconfirmed, and the same error has also been
reported independently on a real Mac (and once under WSL) in
`hklarner/pyboolnet#99 <https://github.com/hklarner/pyboolnet/issues/99>`_,
with no root cause ever identified there either.

Aside from ``steady_states()``, BoolDog is unaffected by the gringo/clasp
issue: model instantiation, Boolean/continuous simulation, and I/O all
work fine once Rosetta 2 is installed. Neither issue is something BoolDog
can fix directly, since both are upstream pyboolnet limitations; see
`NIB-SI/BoolDog#6 <https://github.com/NIB-SI/BoolDog/issues/6>`_ for
details and status.

Alternative installations
=========================

Google Colab
------------

BoolDog can be used via
`Google Colab <https://colab.research.google.com/drive/1uUdGZF0Uad6ckkTe3bwZKrzbKpOWsbHB?usp=sharing>`_


Docker
------

A ``Dockerfile`` is provided at the repository root, building on
``python:3.12-slim`` and installing BoolDog with all optional
dependencies.

.. code-block:: bash

   git clone https://github.com/NIB-SI/BoolDog.git
   cd BoolDog
   docker build -t booldog .
   docker run -it --rm booldog python

On an Apple Silicon (M-series) Mac, pass ``--platform linux/amd64`` to
both commands, so Docker emulates x86_64 rather than picking a native
arm64 base image, for which pyboolnet has no matching binaries at all (see
`Docker's multi-platform build docs <https://docs.docker.com/build/building/multi-platform/>`_):

.. code-block:: bash

   docker build --platform linux/amd64 -t booldog .
   docker run --platform linux/amd64 -it --rm booldog python
