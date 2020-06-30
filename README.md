# squad-reboot

Python based squad algorithm. 

## Usage 

See test-notebook.ipynb. 

```python
import squad_reboot
g = squad_reboot.SquadRegulatoryNetwork(squad_reboot.import_graphml("./examples/Athaliana.graphml"))
g.dynamic_simulation(t_max=30, gamma=1, h=10)
```

---

## Dependencies: 

### Not included

* numpy
* xmltodict
* scipy
* python-igraph
* matplotlib

### Included 

* PyBoolNet
