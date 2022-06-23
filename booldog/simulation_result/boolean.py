import numpy as np
import os
import sys
import math

import importlib

from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx

from pyboolnet.state_transition_graphs import stg2image

import warnings


# path to mpl stylesheet
_style_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           'utils', 'stylesheet.mplstyle')


class BooleanSimulationResult():
    """Small class to contain simulation results"""

    def __init__(self, boolean_graph, stg):
        # TODO copy RN so that edits to original do not affect this object
        self.boolean_graph = boolean_graph

        self.stg = stg


    def plot(self, file=None,
           booldog_style=True,
           plot_nodes=None):
        '''Plot the state transition graph, from optional initial values.

        Parameters
        ----------
        fig : str
            File name of generated figure

        initial_values : int, str, list, dict
            Initial state, see Notes for format

        new_style : bool
            Whether to use booldog style, or PyBoolNet style (default) to
            plot the state transition graph. Requires pygraphviz (If not
            installed, will default to PyBoolNet version.)

        plot_nodes :  None or list of str, optional
            Subset of nodes to plot. If `None`, plot all nodes.
            Only valid if `booldog_style` is `True`.
        '''

        if booldog_style and (importlib.util.find_spec('pygraphviz') is not None):

            if plot_nodes is None:
                plot_nodes = self.boolean_graph.nodes

            rev_index = {}
            for node in self.boolean_graph.nodes:
                rev_index[self.boolean_graph.index[node]] = node

            self.stg.graph['node']['shape'] = 'plaintext'
            colors = {"1":"#80ff8a", "0":"#ff9580"}

            for n in self.stg:
                label = '<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="1" ><TR>'
                for i, x in enumerate(n):
                    node_name = rev_index[i]
                    if node_name in plot_nodes:
                        label += f'<TD BGCOLOR="{colors[x]}">{rev_index[i][:5]} <BR/> {x} </TD>'
                label += """</TR></TABLE>>"""
                self.stg.nodes[n]["label"]=label

            A = nx.drawing.nx_agraph.to_agraph(self.stg)
            A.layout('dot')

            if file:
                A.draw(file)
                print(f"Saved figure to {file}. ")

            return A

        else:
            if booldog_style:
                warnings.warn("pygraphviz not available, defaulting to pyboolnet style.", ImportWarning)

            self.stg.graph["node"]["color"] = "cyan"
            self.stg.graph["node"]["height"] = 0.3
            self.stg.graph["node"]["width"] = 0.45

            if file is not None:
                stg2image(self.stg, file, layout_engine="dot")

            return self.stg.graph