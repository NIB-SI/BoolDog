import numpy as np
import os
import sys
import math

from pathlib import Path

import matplotlib.pyplot as plt


# path to mpl stylesheet
_style_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                        'stylesheet.mplstyle')



class ContinousSimulationResult():
    """Small class to contain simulation results"""

    def __init__(self, regulatory_network, t, y, ode_system,
                 node_events=None, edge_events=None
        ):

        # copy RN so that edits to original do not affect this object
        self.regulatory_network = regulatory_network#.copy()

        self.t = t
        self.y = y

        self.node_events = node_events
        self.edge_events = edge_events

    def export(self, outfile):

        # check if expert path is "writeable" if is not False:
        outfile = Path(outfile)
        file_writeable(outfile)

        with open(outfile, "w") as out:
            # write node list (for order)
            out.write("#nodelist\t" + "\t".join(self.regulatory_network.nodes) + "\n")

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

        print(f"Saved simulation results to {outfile}. ")

    def plot(self, file=None,
             plot_nodes=None,
             title=None,
             figsize=(20, 10)
        ):
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
            vlines += [x['time'] + x['duration'] for x in self.node_events if 'duration' in x]
        if self.edge_events:
            vlines += self.edge_events

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
                legend_labels = self.regulatory_network.index
                self._plot_one_ax(axes[0, 0], self.t, self.y,
                             legend_labels,
                             vlines=vlines,
                             title=title)
            # case (2) and (3)
            else:
                if not isinstance(plot_nodes[0], list):
                    # case (2) --> case (3)
                    plot_nodes=[plot_nodes]

                num_plots = len(plot_nodes)

                if title and isinstance(title, str):
                    main_title = title
                    title = [None]*num_plots
                elif  title and isinstance(title, list):
                    if not (len(title) == len(plot_nodes)):
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
                        yidx.append(self.regulatory_network.index[node])
                        legend_labels.append(node)
                    this_y = self.y[:, yidx]

                    ax = axes[i, 0]
                    self._plot_one_ax(ax, self.t, this_y,
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


        if file is not None:
            plt.savefig(file, bbox_inches = "tight")

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
