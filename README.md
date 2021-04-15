# squad-reboot

Python based squad algorithm.

## Installation

To install:

```bash
git clone https://github.com/NIB-SI/squad-reboot.git
cd squad-reboot
pip install . 
```

To remove:

```bash
pip uninstall squad_reboot
```

### Dependencies:

* numpy
* xmltodict
* scipy
* python-igraph
* matplotlib
* PyBoolNet
* pygraphviz

## Usage

See test-notebook.ipynb.

```python
import squad_reboot
g = squad_reboot.RegulatoryNetwork("./examples/Athaliana.graphml", "graphml")
g.continuous_simulation(t_max=30, gamma=1, h=10)
```

---


