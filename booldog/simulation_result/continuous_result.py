'''Continuos simulation result class'''

from pathlib import Path

import math
import matplotlib.pyplot as plt

from booldog.resources import mpl_style_sheet
from booldog.utils import file_writable, get_pkg_version

import logging

logger = logging.getLogger(__name__)

class ContinuousSimulationResult():
    """Class to contain simulation results"""

    def __init__(self,
                 network,
                 t,
                 y,
                 ode_system,
                 node_events=None,
                 edge_events=None):

        # copy Network so that edits to original do not affect this object
        self.network = network  #.copy()
        self.ode_system = ode_system

        self.t = t
        self.y = y

        self.node_events = node_events
        self.edge_events = edge_events

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
        If you want to use `pandas` to read the file, you can use the following code:
        df = pd.read_csv(outfile, sep="\t")
        '''

        # check if expert path is "writeable" if is not False:
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
        """Plot simulation results.

        Called by :py:func:`continous_simulation`.

        Parameters
        ----------
        plot_nodes :  None or list of str or list of lists of str, optional
            Subset of nodes to plot. If `None`, plot all nodes. If a list of
            lists, each sublist is plotted as a subplot.
        title : None or string or list of string
            If str, main title of the plot. If a list of str, subtitles of
            subplots as defined by `plot_node`. In this case `plot_nodes`
            should be a list of lists, and `title` should be the same length as
            `plot_nodes`.
        figsize : (float, float)
            Width, height in inches.

        Returns
        ----------
        fig :  matplotlib.figure.Figure
        axes : array of matplotlib.axes.Axes

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
