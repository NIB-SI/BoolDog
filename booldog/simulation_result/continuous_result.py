'''Continuous simulation result class.

Contains :class:`ContinuousSimulationResult`, returned by
:py:meth:`~booldog.continuous.semi_quantitative.ContinuousMixin.continuous_simulation`.
'''

from pathlib import Path

import math
import matplotlib.pyplot as plt

from booldog.resources import mpl_style_sheet
from booldog.utils import file_writable, get_pkg_version

import logging

logger = logging.getLogger(__name__)

class ContinuousSimulationResult():
    """Class to contain the result of a continuous (semi-quantitative) ODE
    simulation of a Boolean network: the time-points and solution values,
    the :class:`~booldog.continuous.ode_factory.ODE` system used to generate
    them, and any node/edge perturbation events applied. Returned by
    :py:meth:`~booldog.continuous.semi_quantitative.ContinuousMixin.continuous_simulation`."""

    def __init__(self,
                 network,
                 t,
                 y,
                 ode_system,
                 node_events=None,
                 edge_events=None):
        '''
        Parameters
        ----------
        network : BoolDogModel
            The network the simulation was run on.
        t : ndarray
            1D array of time-points of the simulation.
        y : ndarray
            2D array of simulated values, of shape
            ``(len(t), len(network.nodes))``; column order matches
            `network.index`/`network.node_ids`.
        ode_system : ODE
            The :class:`~booldog.continuous.ode_factory.BooleCubeODE` or
            :class:`~booldog.continuous.ode_factory.SquadODE` instance used
            to generate the simulation.
        node_events : list of dict or None, optional
            Node perturbation events applied during the simulation, each a
            dict with keys ``'time'``, ``'node'``, ``'value'`` and
            (optionally) ``'duration'`` — see
            :py:meth:`~booldog.continuous.semi_quantitative.ContinuousMixin.continuous_simulation`
            for details.
        edge_events : list or None, optional
            Edge perturbation events; not currently implemented upstream
            (see Notes).
        '''

        # copy Network so that edits to original do not affect this object
        self.network = network  #.copy()
        '''BoolDogModel : the network the simulation was run on.'''

        self.ode_system = ode_system
        '''ODE : the :class:`~booldog.continuous.ode_factory.BooleCubeODE` or
        :class:`~booldog.continuous.ode_factory.SquadODE` instance used to
        generate this simulation.'''

        self.t = t
        '''ndarray : 1D array of time-points of the simulation.'''

        self.y = y
        '''ndarray : 2D array of simulated values, shape
        ``(len(t), len(network.nodes))``, column order matching
        `network.index`/`network.node_ids`.'''

        self.node_events = node_events
        '''list of dict or None : node perturbation events applied during the
        simulation (see `__init__` Parameters).'''

        self.edge_events = edge_events
        '''list or None : edge perturbation events; not currently implemented
        upstream.'''

    def export(self, outfile, decimals=5):
        '''
        Export simulation results to a file.

        Parameters
        ----------
        outfile : str or Path
            Path to the output file. The file will be created if it does not
            exist. If the file already exists, it will be overwritten.

        decimals : int
            Number of decimals to round the output values to. Default is 5.

        Notes
        -----
        The output file will contain:

            - nodelist
            - ODE transform
            - ODE parameters
            - node_events
            - timepoints
            - solution/y

        The output will be tab-separated and can be read into a pandas DataFrame.
        If you want to use `pandas` to read the file, you can use the following code::

            df = pd.read_csv(outfile, sep="\\t")

        The "ODE parameters" line(s) are written from
        ``ode_system.param_dict``, which both
        :class:`~booldog.continuous.ode_factory.BooleCubeODE` and
        :class:`~booldog.continuous.ode_factory.SquadODE` provide.
        '''

        # check if export path is "writeable" if is not False:
        outfile = Path(outfile)
        file_writable(outfile)

        with open(outfile, "w", encoding="utf-8") as out:
            # write file origin
            out.write(f"#Semi-quantitative simulation results exported from booldog version {get_pkg_version()}.\n")

            # write model source
            if self.network.modelinfo.source:
                out.write(f"#Model source: {self.network.modelinfo.source}.\n")

            # write node list (for order)
            out.write("#nodelist\t" +
                      "\t".join(self.network.nodes) + "\n")

            # write transform
            out.write("#transform\t" + self.ode_system.transform + "\n")

            # write parameters
            for param_name, param in self.ode_system.param_dict.items():
                out.write("#param\t" + param_name + "\t"+ \
                          "\t".join(param.astype(str)) + "\n")

            # write events
            out.write("#node_events\t" + str(self.node_events) + "\n")

            # write timepoints
            out.write("time\t" + \
                      "\t".join(self.t.round(decimals).astype(str)) + "\n")

            # write solution/y
            for node, array in zip(self.network.nodes, self.y.T):
                out.write(node + "\t" + \
                          "\t".join(array.round(decimals).astype(str)) + "\n")

        logger.info("Saved simulation results to %s." , outfile)

        # time 1 2 3 4 5
        # x x1 x2 x3
        # z z1 z2 z3


    def plot(self, file=None, plot_nodes=None, title=None, figsize=(20, 10)):
        """Plot the simulated time-series for each node.

        Each requested node (or group of nodes, see `plot_nodes`) is plotted
        as a line of relative concentration against time, on axes shared
        across subplots. Vertical dashed lines mark the start (and, if given,
        end) times of any `node_events`/`edge_events` recorded on this
        result.

        Parameters
        ----------
        file : str or Path or None, optional
            If given, the figure is saved to this path (via
            ``plt.savefig(file, bbox_inches="tight")``) instead of being
            shown interactively. Default `None`.
        plot_nodes :  None or list of str or list of lists of str, optional
            Subset of nodes to plot. If `None`, plot all nodes on a single
            axes. If a list of node identifiers, plot only those nodes on a
            single axes. If a list of lists of node identifiers, each sublist
            is plotted on its own subplot (one row per sublist).
        title : None or str or list of str, optional
            If `plot_nodes` is `None`/empty and `title` is a str, it is used
            as the (single) axes title. If `plot_nodes` is given (flat or
            nested list) and `title` is a str, it is used as a figure-level
            title instead (see Notes). If a list of str, used as the
            per-subplot subtitles as defined by `plot_nodes`; `plot_nodes`
            should then be a list of lists, and `title` should be the same
            length as `plot_nodes` (otherwise a warning is logged and
            subtitles are omitted).
        figsize : (float, float), optional
            Width, height in inches, passed to ``plt.subplots``. Default
            ``(20, 10)``.

        Returns
        -------
        fig : matplotlib.figure.Figure
        axes : ndarray of matplotlib.axes.Axes
            Array of axes, one row per (sub)plot.

        Notes
        -----
        If `plot_nodes` is `None`/empty, a single axes is drawn and a string
        `title` (if given) is used directly as that axes' title. Otherwise
        (`plot_nodes` given, flat or nested), a string `title` is instead used
        as a figure-level ``suptitle``, and per-axes titles are only set from
        a list `title` matching the number of (sub)plots.
        """

        # collect vertical lines at events
        vlines = []
        if self.node_events:
            vlines += [x['time'] for x in self.node_events]
            vlines += [
                x['time'] + x['duration'] for x in self.node_events
                if 'duration' in x
            ]
        if self.edge_events:
            vlines += self.edge_events

        with plt.style.context(mpl_style_sheet):
            # 3 cases to plot:
            # (1) no node list
            # (2) one node list
            # (3) multiple node lists

            main_title = False

            # case (1)
            if not plot_nodes:
                fig, axes = plt.subplots(ncols=1,
                                         nrows=1,
                                         squeeze=False,
                                         figsize=figsize,
                                         constrained_layout=True)
                legend_labels = [self.network.nodes[node_id].name for node_id in self.network.index]
                self._plot_one_ax(axes[0, 0],
                                  self.t,
                                  self.y,
                                  legend_labels,
                                  vlines=vlines,
                                  title=title)
            # case (2) and (3)
            else:
                if not isinstance(plot_nodes[0], list):
                    # case (2) --> case (3)
                    plot_nodes = [plot_nodes]

                num_plots = len(plot_nodes)

                if title and isinstance(title, str):
                    main_title = title
                    title = [None] * num_plots
                elif title and isinstance(title, list):
                    if not (len(title) == len(plot_nodes)):
                        logger.warning(
                              'Number of (sub)titles is not equal to the '\
                              'number of (sub)plots. Either pass the correct '\
                              'number of subtitles as a list, or a single '\
                              'main title as a string. \n'\
                              ' %s (title) != %s (subplots)', len(title), num_plots)
                        title = [None] * num_plots

                else:
                    title = [None] * len(plot_nodes)

                fig, axes = plt.subplots(ncols=1,
                                         nrows=num_plots,
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
                    for node_id in node_list:
                        yidx.append(self.network.index[node_id])
                        legend_labels.append(self.network.nodes[node_id].name)
                    this_y = self.y[:, yidx]

                    ax = axes[i, 0]
                    self._plot_one_ax(ax,
                                      self.t,
                                      this_y,
                                      legend_labels,
                                      vlines=vlines,
                                      title=subtitle)

            fig.supxlabel('Time', fontsize=18, fontweight='bold')
            fig.supylabel("Relative concentration",
                          x=-0.02,
                          fontsize=18,
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

        if file is not None:
            plt.savefig(file, bbox_inches="tight")
        else:
            plt.show()

        return fig, axes

    def _plot_one_ax(self, ax, x, y, legend_labels, vlines=None, title=None):
        '''Helper method that draws one subplot of :py:meth:`plot`: plots
        `y` against `x` as lines with a legend, fixed y-limits of (0, 1),
        offset spines, and optional vertical event markers.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axes to plot on.
        x : ndarray
            1D array of time-points (x-values), shared by all lines.
        y : ndarray
            2D array of values to plot, one column per line/node (passed
            directly to ``ax.plot(x, y)``).
        legend_labels : list of str
            Legend label for each column of `y`, in the same order.
        vlines : list of float or None, optional
            x-positions at which to draw dashed grey vertical lines (e.g.
            event times). Default `None` (no lines).
        title : str or None, optional
            If given, set as this axes' title.

        Returns
        -------
        None
            The axes `ax` is modified in place.
        '''

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
                      ymin=ymin,
                      ymax=ymax,
                      colors='gray',
                      ls='--',
                      alpha=0.5)

        lines = ax.plot(x, y)
        legend = ax.legend(lines,
                           legend_labels,
                           bbox_to_anchor=(1.01, 0.5),
                           ncol=math.ceil(len(legend_labels) / 20))
        legend.set_in_layout(False)

        if title:
            ax.set_title(title)
