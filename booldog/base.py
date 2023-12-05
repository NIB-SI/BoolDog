import numpy as np
import os
import sys
import math

from pathlib import Path

# plotting library #TODO optional import?
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
from .ode import ODE_factory

# path to mpl stylesheet
_style_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           'utils', 'stylesheet.mplstyle')

class RegulatoryNetwork(BooleanGraph):
    '''A class to represent a regulatory network.

    '''

    def __init__(self, graph, **kwargs):
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
        super().__init__(graph,  **kwargs)

    def transform_bool_to_continuous(self, transform, **kwargs):
        '''Generate an ODE from RegulatoryNetwork/Boolean graph.

        Parameters
        ----------
        transform : str
            One of accepted transforms. See `booldog.ode.transforms` for
            accepted options.

        Other Parameters
        ----------
        **kwargs
            Additional keyword arguments passed to :py:func:`booldog.ODE`.
        '''
        return ODE_factory(self, transform, **kwargs)

    def continuous_simulation(self,
        node_events={},
        edge_events={},
        t_min=0, t_max=30,
        initial_state=0,
        plot=True,
        export=False,
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
        plot : bool, optional
            Whether to plot the simultion results
        export : bool or path object or string, optional
            False, or a file path to save simulation results.
            Exports values to 5 decimal points.
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

        # check if expert path is "writable" if is not False:
        if export:
            export = Path(export)
            file_writeable(export)

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

        # add a last event to get to end of simulation
        # i.o.t. avoid using `while True`
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

        if plot:
            try:
                self.plot_simulation(combined_t, combined_y,
                                node_events=sorted(events_d.keys())[:-1],
                                **kwargs)
            except Exception as e:
                print("Error in plotting")
                print(e)

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

    def plot_simulation(self,
        t, y,
        node_events=None,
        edge_events=None,
        plot_nodes=None,
        title=None,
        figsize=(20, 10),
        **kwargs):
        """Plot simulation results.

        Called by :py:func:`continous_simulation`.

        Parameters
        ----------
        t : ndarray, shape (n_time_points,)
            Time-points.
        y : ndarray, shape (n_time_points, n_nodes)
            Values of the solution at t.
        plot_nodes :  None or list of str or list of lists of str, optional
            Subset of nodes to plot. If `None`, plot all nodes. If a list of
            lists, each sublist is plotted as a subplot.
        node_events : None or list of dict, optional
            List of node events with a dictionary defining each event.
            See :ref:`Notes <tagnotesne>` for description of event definitions.
        edge_events :  None or list of dict, optional
            Disrupt connections #TODO not implemented
        title : None or string or list of string
            If str, main title of the plot. If a list of str, subtitles of
            subplots as defined by `plot_node`. In this case `plot_nodes`
            should be a list of lists, and `title` should be the same length as
            `plot_nodes`.
        figsize : (float, float)
            Width, height in inches.

        Other Parameters
        ----------
        **kwargs
            ignored

        Returns
        ----------
        fig :  matplotlib.figure.Figure
        axes : array of matplotlib.axes.Axes

        """

        # collect vertical lines at events
        vlines = []
        if node_events:
            vlines += node_events
        if edge_events:
            vlines += edge_events

        with plt.style.context(_style_path):
            # 3 cases to plot:
            # (1) no node list
            # (2) one node list
            # (3) multiple node lists

            main_title = False

            # case (1)
            if not plot_nodes:
                fig, axes = plt.subplots(ncols=1, nrows=1,
                                         squeeze=False,
                                         figsize=figsize,
                                         constrained_layout=True)
                legend_labels = self.index
                self._plot_one_ax(axes[0, 0], t, y,
                             legend_labels,
                             vlines=vlines,
                             title=title)
            else:
                if not isinstance(plot_nodes[0], list):
                    # case (2) --> case (3)
                    plot_nodes=[plot_nodes]

                num_plots = len(plot_nodes)

                if title and isinstance(title, str):
                    main_title = title
                    title = [None]*num_plots
                elif title and isinstance(title, list):
                    if len(title) != len(plot_nodes):
                        print(f'WARNING: '\
                                  f'Number of (sub)titles is not equal to the '\
                                  f'number of (sub)plots. Either pass the correct '\
                                  f'number of subtitles as a list, or a single '\
                                  f'main title as a string. \n'\
                                  f'{len(title)} (title) != {num_plots} (subplots)')
                        title = [None]*num_plots

                else:
                    title = [None]*len(plot_nodes)

                fig, axes = plt.subplots(ncols=1, nrows=num_plots,
                                         sharex=True,
                                         sharey=True,
                                         squeeze=False,
                                         figsize=figsize,
                                         constrained_layout=True)

                for i in range(num_plots):
                    node_list = plot_nodes[i]
                    subtitle = title[i]

                    yidx = []
                    legend_labels = []
                    for node in node_list:
                        yidx.append(self.index[node])
                        legend_labels.append(node)
                    this_y = y[:, yidx]

                    ax = axes[i, 0]
                    self._plot_one_ax(ax, t, this_y,
                                 legend_labels,
                                 vlines=vlines,
                                 title=subtitle)

            fig.supxlabel('Time', fontsize=18, fontweight='bold')
            fig.supylabel("Relative concentration", x=-0.02, fontsize=18,
                           fontweight='bold')
            if main_title:
                fig.suptitle(main_title)

            # issues with keeping legend in figure bounds
            # when using constrained_layout see
            # https://matplotlib.org/stable/tutorials/intermediate/constrainedlayout_guide.html#legends #noqa
            fig.canvas.draw()
            for ax in axes.flatten():
                legend = ax.get_legend()
                legend.set_in_layout(True)
            fig.set_constrained_layout(False)
        return fig, axes

    def _plot_one_ax(self,
        ax,
        x, y,
        legend_labels,
        vlines=None,
        title=None):
        '''Helper func that plots on a subplot axes'''

        ymin, ymax = 0, 1
        xmin, xmax = 0, max(x)
        ax.set_ylim((ymin, ymax))
        ax.set_xlim((xmin, xmax))

        # offset spines and limit
        ax.spines['left'].set_position(('outward', 5))
        ax.spines['left'].set_bounds((ymin, ymax))
        ax.spines['bottom'].set_position(('outward', 5))
        ax.spines['bottom'].set_bounds(0, int(xmax))

        if vlines:
            ax.vlines(x=vlines,
                   ymin=ymin, ymax=ymax,
                   colors='gray', ls='--', alpha=0.5)

        lines = ax.plot(x, y)
        legend = ax.legend(lines, legend_labels,
                  bbox_to_anchor=(1.01, 0.5),
                  ncol=math.ceil(len(legend_labels)/20)
                  )
        legend.set_in_layout(False)

        if title:
            ax.set_title(title)