import numpy as np
import os
import sys

from pathlib import Path

import matplotlib.pyplot as plt

from collections import defaultdict
from scipy.integrate import solve_ivp

# hide runtime warnings (divide by zero, multiply by inf)
# TODO This is a bad idea...
import warnings
warnings.filterwarnings("ignore", message="divide by zero encountered in true_divide")
warnings.filterwarnings("ignore", message="invalid value encountered in multiply")

from .utils import *
from .bool_graph import BooleanGraph
from .ode import *

class RegulatoryNetwork:
    '''
    A class to represent a regulatory network.

    Attributes
    ----------
    n : int
        The number of nodes/variables in the network.
    boolean_graph : BooleanGraph
        The associated Boolean graph representation of the network.
    is_interaction : Bool
        If the graph is an interaction network, and not a Boolean graph.

    Methods
    ----------
    boolean_simulation

    transform_bool_to_continuous

    continuous_simulation


    '''

    def __init__(self, graph, **kwargs):
        '''
        Initialise a regulatory network #TODO from a Boolean graph.

        Parameters
        ----------
        graph : squad_reboot.BooleanGraph, dict or str
            If not a BooleanGraph object, then correct type for initalising
            a squad_rebbot.BooleanGraph

        Other Parameters
        ----------
        **kwargs
            Additional arguments and keyword arguments passed to squad_reboot.BooleanGraph.
            For description of the arguments see help(squad_reboot.BooleanGraph).

        '''
        if isinstance(graph, BooleanGraph):
            self.boolean_graph = graph

        else:
            self.boolean_graph = BooleanGraph(graph,  **kwargs) # TODO add format translations?

        self.n = self.boolean_graph.n




    def boolean_simulation(self, initial_state, num_transitions):
        #TODO
        pass

    def transform_bool_to_continuous(self, transform, **kwargs):
        transform = transform.lower()
        if not transform in ode_classes.keys():
            raise ValueError("'transform' argument must be one of %s"%str(list(ode_classes.keys())))

        ode_system = ode_classes[transform](self.boolean_graph, transform, **kwargs)
        return ode_system

    def continuous_simulation(self,
       node_events={},
       edge_events={},
       t_min=0, t_max=30,
       initial_state=0,
       plot=True,
       plot_nodes=[],
       export=False,
       ode_system=None,
       **kwargs):
        '''
        Run continuous qualitative simulation.

        Parameters
        ----------
        node_events : dict, optional
            Dictionary of node events with a key, value pair for each time point at
            which one or more events occur:
                key - time at which the event occurs
                value - dictionary with affected node identifier as key and
                        values defining the event ()
            Example - at timepoint 10, node X is set to 0.25 for 5 time-steps.
                      and at timepoint 12, node Y and Z are set to 1 for zero timesteps:
                            node_events = {
                                10:{"X":{"value":.25, "duration":5, "perturbation":"set"}},
                                20:{"Y":{"value":1, "perturbation":"set"},
                                    "Z":{"value":1, "perturbation":"set"}
                                   }
                                }


        edge_events : dict, optional
            Disrupt connections

        t_min : float, optional

        t_max : float, optional

        initial_state : dict, optional

        plot : bool, optional

        plot_nodes : list, optional

        export : bool, path object or string
            False, or path to save simulation results

        ode_system: None or ODE class


        Other Parameters
        ----------
        **kwargs
            If ode_system is None, additional keyword arguments are
            passed to transform_bool_to_continuous.

            For description of these arguments see help(squad_reboot.ODE).

        Returns
        ----------
        t : numpy array
            Time-points.
        y : numpy array
            Value of the solution at time-points.

        '''

        # check if expert path is "writable" if is not False:
        if export:
            export = Path(export)
            file_writeable(export)

        if ode_system is None:
            # fetch the system of eqn
            ode_system = self.transform_bool_to_continuous(**kwargs)

        all_results = []

        initial_state_array = np.zeros(self.n)
        if not (initial_state == 0):
            for node, value in initial_state.items():
                initial_state_array[self.boolean_graph.index[node]] = value

        # 1. copy node_events to new dict
        # 2. set a duration for all node_events
        # 3. add end of each perturbation as an event (so that dx/dt of the node can be reset)
        events_d = defaultdict(dict)
        for event_t  in node_events.keys():
            for node, perturb in node_events[event_t].items():
                events_d[event_t][node] = perturb
                events_d[event_t][node]["reset"] = False
                if not ("duration" in perturb):
                    perturb_duration = events_d[event_t][node]["duration"] = 0
                else:
                    perturb_duration = perturb["duration"]

                if perturb_duration > 0:
                    perturb_end = min(t_max, int(event_t + perturb_duration))
                    events_d[perturb_end] = {node:{"reset":True}}

        # add a last event to get to end of simulation
        # i.o.t. avoid using `while True`
        # # https://stackoverflow.com/a/50703835/4996681
        events_d[t_max] = {}

        off_nodes = set()
        for event_t  in sorted(events_d.keys()):
            print("Status: Start")

            result = solve_ivp(fun=ode_system.dxdt, t_span=(t_min, t_max), y0=initial_state_array,
                                         events=ode_system.event_function, args=(event_t,),
                                         dense_output=True, max_step=0.01)
            all_results.append(result)
            current_t = result.t[-1]

            if (result.status == 0) or (current_t == t_max):
                print("Status: End")
                break
            else:
                s = "Status: Event at {t:3d}:\n".format(t=int(current_t))

                t_min = current_t
                initial_state_array = result.y[:, -1].copy()
                for node, perturb in events_d[event_t].items():

                    if perturb["reset"]:
                        # ending a perturbation (put back dx/dt of this node)
                        s += " "*10 + "{node} -> released\n".format(node=node)
                        off_nodes.remove(self.boolean_graph.index[node])
                        ode_system.update(off=off_nodes)

                    else:
                        # starting a perturbation
                        initial_state_array[self.boolean_graph.index[node]] = perturb['perturbation']
                        s += " "*10 + "{node} -> {perturb_value:.2f} (duration {duration:2d})\n".format(\
                                                           node=node,
                                                           perturb_value=perturb['perturbation'],
                                                           duration=perturb["duration"])

                        if perturb['duration'] > 0:
                            # dt/dt of this node should be 0 for "duration"
                            off_nodes.add(self.boolean_graph.index[node])
                            ode_system.update(off=off_nodes)
                print(s)

        combined_t = np.concatenate([result.t for result in all_results])
        combined_y = np.concatenate([result.y for result in all_results], axis=1).T

        if plot:

            if len(plot_nodes) == 0:
                y_plot = combined_y
                legend_labels = self.boolean_graph.index
            else:
                yidx = []
                legend_labels = []
                for node in plot_nodes:
                    yidx.append(self.boolean_graph.index[node])
                    legend_labels.append(node)
                y_plot = combined_y[:, yidx]


            lines = plt.plot(combined_t, y_plot, '-')
            plt.legend(lines, legend_labels)
            plt.ylim([-0.01, 1.01])

        if export:
            pass
            #TODO


        return combined_t, combined_y