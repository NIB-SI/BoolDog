import numpy as np
import os
import sys

from pathlib import Path

import matplotlib.pyplot as plt

from collections import defaultdict
from scipy.integrate import solve_ivp

# hide runtime warnings (divide by zero, multiply by inf)
# occur in squad ode init
# TODO This is a bad idea...
import warnings
warnings.filterwarnings("ignore",
    message="divide by zero encountered in true_divide")
warnings.filterwarnings("ignore",
    message="invalid value encountered in multiply")

from .utils.utils import *
from .bool_graph import BooleanGraph
from .ode import ODE, ode_classes

class RegulatoryNetwork(BooleanGraph):
    '''
    A class to represent a regulatory network.
    Inherits from BooleanGraph.

    '''

    def __init__(self, graph, **kwargs):
        '''
        Initialise a regulatory network from a Boolean graph.

        Parameters
        ----------
        graph : booldog.BooleanGraph, dict or str
            If not a BooleanGraph object, then correct input for initializing
            a booldog.BooleanGraph.

        Other Parameters
        ----------
        **kwargs
            Additional arguments and keyword arguments passed to
            booldog.BooleanGraph. For description of these arguments see
            help(booldog.ODE).
        '''
        super().__init__(graph,  **kwargs)

    def transform_bool_to_continuous(self, transform, **kwargs):
        '''
        Generate an ODE from RegulatoryNetwork/Boolean graph.

        Parameters
        ----------
        transform : str
            One of accepted transforms. See `booldog.ode_transforms` for
            options.

        Other Parameters
        ----------
        **kwargs
            Additional arguments and keyword arguments passed to
            booldog.ODE. For description of these arguments see
            help(booldog.ODE).
        '''

        transform = transform.lower()
        if not transform in ode_classes.keys():
            raise ValueError(f"transform' argument must be one of"\
                             f"{list(ode_classes.keys())}")

        ode_system = ODE(self, transform,**kwargs)
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
        '''Run continuous semi-qualitative simulation.

        Parameters
        ----------
        node_events : list of dict, optional
            List of node events with a dictionary defining each event.See notes for description of event definitions.

        edge_events : dict, optional
            Disrupt connections #TODO not implemented

        t_min : float, optional

        t_max : float, optional

        initial_state : dict, optional

        plot : bool, optional

        plot_nodes : list, optional
            If list is non-empty, only plot simulation results of listed
            nodes. Otherwise plot all nodes.

        export : bool, path object or string
            False, or path to save simulation results.
            Exports values to 5 decimal points.
            #TODO describe export format.

        ode_system: ODE class object, optional
            If not set, is created


        Other Parameters
        ----------
        **kwargs
            If ode_system is None, additional keyword arguments are
            passed to transform_bool_to_continuous.

            For description of these arguments see help(booldog.ODE).

        Returns
        ----------
        t : numpy array
            Time-points.
        y : numpy array
            Value of the solution at time-points.

        Notes
        -----
        Format of the `node_events` parameter:
        The node events are passed as a list of dictionaries defining each event.
        Dictionary keys are:
            time: time at which the event occurs
            node: node which is perturbed
            value:  value at which the node is set
            duration: (optional) duration for which the node is fixed
                    if longer than 1, (i.e. not a point perturbation)

        Example - at timepoint 10, node X is set to 0.25 for 5 time-steps.
          and at timepoint 12, node Y and Z are set to 1 for 1
          timesteps::

            node_events = [
                {'time':10, 'node':'X', 'value':.25, 'duration':5},
                {'time':12, 'node':'Y', 'value':1},
                {'time':12, 'node':'X', 'value':1}
            ]

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
                initial_state_array[self.index[node]] = value

        # 1. copy node_events to new dict
        # 2. set a duration for all node_events
        # 3. add end of each perturbation as an event
        #    (so that dx/dt of the node can be reset)
        events_d = defaultdict(lambda: defaultdict(dict))

        for node_event  in node_events:

            this_time_and_node = {
                    "reset":False,
                    "duration":0,
                    "value":node_event["value"]
            }
            if ( "duration" in node_event ) and ( node_event["duration"] > 0 ):
                this_time_and_node["duration"] = node_event["duration"]

                perturb_end = int(node_event['time'] + node_event["duration"])
                events_d[perturb_end][node_event['node']]["reset"] = True

            events_d[node_event['time']][node_event['node']] = this_time_and_node

        # add a last event to get to end of simulation
        # i.o.t. avoid using `while True`
        # # https://stackoverflow.com/a/50703835/4996681
        events_d[t_max] = {}

        for time  in sorted(events_d.keys()):
            print(time)
            for node, event in events_d[time].items():
                print(f"    {node} {event}")


        print("initial_state: ", initial_state)

        off_nodes = set()
        print("Status: Start")
        for event_t  in sorted(events_d.keys()):

            result = solve_ivp(fun=ode_system.dxdt,
                               t_span=(t_min, t_max),
                               y0=initial_state_array,
                               events=ode_system.event_function,
                               args=(event_t,),
                               max_step=0.01,
                               )

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
                        off_nodes.remove(self.index[node])
                        ode_system.update(off_nodes=off_nodes)

                    else:
                        # starting a perturbation
                        initial_state_array[self.index[node]] =\
                                             perturb['value']

                        s += f"{node:>10} -> {perturb['value']:.2f} "\
                             f"(duration {perturb['duration']:2d})\n"

                        if perturb['duration'] > 0:
                            # dt/dt of this node should be 0 for "duration"
                            off_nodes.add(self.index[node])
                            ode_system.update(off_nodes=off_nodes)
                print(s)

        combined_t = np.concatenate([result.t for result in all_results])
        combined_y = np.concatenate([result.y for result in all_results],
                                    axis=1).T

        if plot:
            ymin, ymax = -0.05, 1.1

            # add vertical lines at events
            plt.vlines(x=sorted(events_d.keys())[:-1],
                       ymin=ymin, ymax=ymax,
                       colors='gray', ls='--', alpha=0.5)
            plt.vlines(x=edge_events.keys(),
                       ymin=ymin, ymax=ymax,
                       colors='azure', ls='--', alpha=0.5)


            if len(plot_nodes) == 0:
                y_plot = combined_y
                legend_labels = self.index
            else:
                yidx = []
                legend_labels = []
                for node in plot_nodes:
                    yidx.append(self.index[node])
                    legend_labels.append(node)
                y_plot = combined_y[:, yidx]


            lines = plt.plot(combined_t, y_plot, '-')
            plt.legend(lines, legend_labels)

            plt.ylim([ymin, ymax])

        if export:
            with open(export, "w") as out:
                # write node list (for order)
                out.write("#nodelist\t" + "\t".join(self.nodes) + "\n")

                # write transform
                out.write("#transform\t" + ode_system.transform + "\n")

                # write parameters
                for param_name, param in ode_system.param_dict.items():
                    out.write("#param\t" + param_name + "\t"+ \
                              "\t".join(param.astype(str)) + "\n")

                # write events
                out.write("#node_events\t" + str(node_events) + "\n")

                # write timepoints
                out.write("time\t" + \
                          "\t".join(combined_t.round(5).astype(str)) + "\n")

                # write solution/y
                for node, array in zip(self.nodes, combined_y.T):
                    out.write(node + "\t" + \
                              "\t".join(array.round(5).astype(str)) + "\n")

            print(f"Saved simulation results to {export}. ")

        return combined_t, combined_y