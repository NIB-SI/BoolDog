import numpy as np
import os
import sys
import math

from pathlib import Path

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

from booldog.utils.utils import *
from booldog.bool_graph import BooleanGraph
from booldog.ode import ODE_factory

from booldog.simulation_result import ContinousSimulationResult

from booldog.io.write import WriteMixin

class RegulatoryNetwork(BooleanGraph, WriteMixin):
    '''A class to represent a regulatory network.

    '''

    def __init__(self, primes, nodes, index, n):
        '''Initialise a regulatory network from a Boolean graph.

        Parameters
        ----------
        graph : `booldog.BooleanGraph` or dict or str
            If not a :py:class:`booldog.BooleanGraph` instance, then correct
            input for initializing a :py:class:`booldog.BooleanGraph`.

        Other Parameters
        ----------
        **kwargs
            In the case `graph` is not a BooleanGraph instance, additional
            keyword arguments passed to :py:class:`booldog.BooleanGraph`.
        '''
        super().__init__(primes, nodes, index, n)

    def transform_bool_to_continuous(self, transform, **kwargs):
        '''Generate an ODE from RegulatoryNetwork/Boolean graph.

        Note that the BooleanGraph object is kept in memory as the primes of the Boolean network. This means that importing the graph may take a while, depending on the size of the network.



        Parameters
        ----------
        transform : str
            One of accepted transforms. See `booldog.ode.transforms` for
            accepted options.

        Other Parameters
        ----------
        **kwargs
            Additional keyword arguments passed to :py:func:`booldog.ODE`.


        Returns
        ----------


        '''
        ode_system = ODE_factory(self, transform, **kwargs)
        return ode_system

    def continuous_simulation(self,
        node_events={},
        edge_events={},
        t_min=0, t_max=30,
        initial_state=0,
        ode_system=None,
        **kwargs):
        '''Run continuous semi-qualitative simulation.

        Parameters
        ----------
        node_events :  None or list of dict, optional
            List of node events with a dictionary defining each event.
            See :ref:`Notes <tagnotesne>` for description of event definitions.
        edge_events :  None or list of dict, optional
            Disrupt connections #TODO not implemented
        t_min, t_max : float, optional
            Interval of integration, simulation starts with `t=t_min` and
            integrates until it reaches `t=t_max`.
        initial_state : float or int or array or dict, optional
            Initial state of nodes. See :ref:`Notes <tagnotesis>` for
            description of format.
        ode_system: None or :py:func:`booldog.ODE`, optional
            If none, the ODE is created with
            :py:func:`transform_bool_to_continuous`

        Other Parameters
        ----------------
        **kwargs
            If ode_system is None, additional keyword arguments are
            passed to :py:func:`transform_bool_to_continuous`.
            For description of these arguments see :py:func:`booldog.ODE`.

            If `plot=True` , additional keyword arguments are
            passed to :py:func:`plot_simulation`.

        Returns
        -------
        r : object #TODO
        t : ndarray, shape (n_time_points,)
            Time-points.
        y : ndarray, shape (n_time_points, n_nodes)
            Values of the solution at t.

        Notes
        -----

        .. _tagnotesne:

        Format of the `node_events` parameter
            The node events are passed as a list of dictionaries defining each
            event. Dictionary keys are:

            - `time`: time at which the event occurs
            - `node`: name of node which is perturbed
            - `value`:  value to which the node is set
            - `duration`: (optional) duration for which the node is fixed \
            if longer than 0, (i.e. not a point perturbation)

            Example - at timepoint 10, node X is set to 0.25 for 5 time-steps.
            and at timepoint 12, node Y and Z are set to 1 for 0
            timesteps::

                node_events = [
                    {'time':10, 'node':'X', 'value':.25, 'duration':5},
                    {'time':12, 'node':'Y', 'value':1},
                    {'time':12, 'node':'X', 'value':1}
                ]

        .. _tagnotesis:

        Format of the `initial_state` parameter
            If the initial state is an int or float, the value is assigned for
            all variables. Otherwise the parameter argument should be a dict
            with keys as node names and values for their initial state. In
            this case, if the initial state is not defined for all nodes, a
            `default` key with the default value should also be present in the
            dict.

        Todo
        ----

        - describe export format.

        '''

        if ode_system is None:
            # fetch the system of eqn
            ode_system = self.transform_bool_to_continuous(**kwargs)

        all_results = []

        initial_state_array = parameter_to_array(initial_state, self.index)

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

        # add a last event to get to end of simulation i.o.t. avoid using `while True`
        # # https://stackoverflow.com/a/50703835/4996681
        events_d[t_max] = {}
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
                        initial_state_array[self.index[node]] = perturb['value']

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

        result = ContinousSimulationResult(self, combined_t, combined_y, ode_system, node_events=node_events, edge_events=edge_events)

        return result


