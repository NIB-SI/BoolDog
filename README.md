# BoolDoG <img src="docs/figures/logo.png" raw=true alt="BoolDoG icon"  width="80" align="right" >

A Python package for analyses of Boolean and semi-qualitative Boolean networks.

## Documentation

[nib-si.github.io/BoolDoG](https://nib-si.github.io/BoolDoG)

## Installation

To install:

```bash
    git clone https://github.com/NIB-SI/BoolDoG.git
    cd BoolDoG
    pip install .
```

To remove:

```bash
pip uninstall booldog
```

### Dependencies:

* numpy
* xmltodict
* scipy
* python-igraph
* matplotlib
* pygraphviz (optional)
* PyBoolNet
* networkx

## Usage

See test-notebook.ipynb.

```python
import booldog
g = booldog.RegulatoryNetwork("./examples/Athaliana.graphml", "graphml")
g.continuous_simulation(t_max=30, gamma=1, h=10)
```

---


