# squad-reboot

Python based squad algorithm. 

## Installation

To install:

```bash
git clone --recurse-submodules https://github.com/NIB-SI/squad-reboot.git
cd squad-reboot
pip install -r requirements.txt
```
To remove: 

```bash
pip uninstall squad_reboot
pip uninstall PyBoolNet
```
### Dependencies: 

#### Not included

* numpy
* xmltodict
* scipy
* python-igraph
* matplotlib

#### Included 

* PyBoolNet

## Usage 

See test-notebook.ipynb. 

```python
import squad_reboot
g = squad_reboot.SquadRegulatoryNetwork(squad_reboot.import_graphml("./examples/Athaliana.graphml"))
g.dynamic_simulation(t_max=30, gamma=1, h=10)
```

---


