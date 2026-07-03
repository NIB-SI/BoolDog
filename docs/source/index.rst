.. BoolDoG documentation master file, created by
   sphinx-quickstart on Wed Apr 14 13:42:00 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

========
BoolDog
========

A Python package for integrated Boolean and semi-quantitative network modelling.

| **Homepage:** https://github.com/NIB-SI/BoolDog
| **Documentation:** https://nib-si.github.io/BoolDog
| **PyPI:** https://pypi.org/project/BoolDog
| **License:** GPL-3.0

Overview
--------

BoolDog provides an integrated, Python-native platform for logic-based modelling, requiring minimal prior knowledge of kinetic or mechanistic parameters. Starting from a regulatory or Boolean network, BoolDog supports the full modelling workflow: from network import and visualisation, through Boolean simulation and attractor analysis, to transformation into continuous ordinary differential equation (ODE) systems and semi-quantitative simulation.

Key features include:

- **Model access:** direct download of model and model metadata from  `BioModels <https://www.biomodels.org/>`_ repository
- **Model import and export:** supports BoolNet, SBML-qual, TabularQual, GraphML, SIF, and native Python formats
- **Network-to-Boolean conversion:** automatic transformation of regulatory networks into Boolean models
- **Visualisation:** built in model visualisation in Cytoscape, further visualisation supported through interoperability with NetworkX and igraph
- **Model modification:** add, remove, or update nodes and rules inplace
- **Boolean analysis** synchronous simulation, state transition graphs, attractor and steady-state analysis via `PyBoolNet <https://github.com/hklarner/pyboolnet>`_
- **Boolean dynamics animation:** generate animated GIFs of Boolean simulations using Cytoscape as a layout engine, allowing custom network layouts and visual styles
- **Semi-quantitative ODE transformation:** conversion of Boolean logic into continuous ODE systems using the SQUAD and normalised Hill cube (ODEfy) schemes
- **Event-based continuous simulation:** simulation of perturbations such as node knockouts or forced   activations at defined time points, with visualisation via matplotlib

.. toctree::
   :maxdepth: 1 
   :caption: Manual and guides

   installation
   example
   formats
   misc

.. toctree::
   :caption: API
   :maxdepth: 2

   api/booldog


Citation
=======

If you use BoolDog in your research, please cite:

   Bleker et al. (2026).
   *BoolDog: integrated Boolean and semi-quantitative network modelling in Python*.
   bioRix.
   https://doi.org/TTTTTT

BoolDog implements ODE transformation schemes from:

- SQUAD:
  Di Cara, A., Garg, A., De Micheli, G., Xenarios, I., & Mendoza, L. (2007).
  *Dynamic simulation of regulatory networks using SQUAD*.
  BMC Bioinformatics, 8(1), 462.
  https://doi.org/10.1186/1471-2105-8-462
- ODEfy:
  Krumsiek, J., Pölsterl, S., Wittmann, D. M., & Theis, F. J. (2010)
  *Odefy -- From discrete to continuous models*.
  BMC Bioinformatics, 11(1), 233.
  https://doi.org/10.1186/1471-2105-11-233

BoolDog relies on PyBoolNet for Boolean analyses:

- Klarner, H., Streck, A., & Siebert, H. (2017).
  *PyBoolNet: a python package for the generation, analysis and visualization of boolean networks*.
  Bioinformatics, 33(5), 770-772.
  https://doi.org/10.1093/bioinformatics/btw682


Indices
=======

* :ref:`genindex`
* :ref:`modindex`   

