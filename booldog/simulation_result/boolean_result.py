'''Boolean simulation results
'''

import contextlib
import importlib
import logging
import tempfile
import time
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import ListedColormap
import networkx as nx
from pyboolnet.state_transition_graphs import stg2image

from booldog.utils import file_writable

# path to style files
from booldog.resources import cytoscape_style_xml

try:
    import py4cytoscape as p4c
    from booldog.utils import cytoscape_utils
    _CYTOSCAPE_AVAILABLE = True
except ImportError as e:
    _CYTOSCAPE_AVAILABLE = False

try:
    from PIL import Image
    _PILLOW_AVALIBLE = True
except ImportError as e:
    _PILLOW_AVALIBLE = False

logger = logging.getLogger(__name__)


DEFAULT_COLOURS = {
    1: "#b2df8a", # a light green,
    0: "#6f6f6f" # a medium grey
}

class BooleanStateSpace():
    """Class to contain a subspace of a Boolean Network."""

    def __init__(self, network, state_space):
        """
        Initialize the BooleanState object.

        Parameters
        ----------
        network : object
            The Boolean network object.
        state_space : list of dict or list of str
            Either: A list of dictionaries representing the state space of the network.
            Each dictionary contains node names as keys and their states (0 to n) as values.
            Or: A list of strings representing the state space of the network,
            where each string is a binary representation of the state
            (e.g. "1010" for a network with 4 nodes).
        """
        self.network = network

        if not isinstance(state_space, list):
            raise ValueError("state_space should be a list.")

        # cannot be empty
        if len(state_space) == 0:
            raise ValueError("state_space cannot be empty.")

        # if list of str, cast to list of dict
        if isinstance(state_space[0], str):
            casted_state_space = []
            for state in state_space:
                if len(state) != len(network.node_ids):
                    raise ValueError(
                        f"State '{state}' has length {len(state)}, but network has {len(network.node_ids)} nodes.")
                if not all(c in "01" for c in state):
                    raise ValueError(
                        f"State '{state}' contains characters other than '0' and '1'.")
                casted_state_space.append(
                    dict(zip(network.node_ids, map(int, state))))
            state_space = casted_state_space

        # if dict, make sure keys are valid node names and values are 0 or 1
        elif isinstance(state_space[0], dict):
            for state in state_space:
                if not all(node in network.node_ids for node in state.keys()):
                    raise ValueError(
                        f"State '{state}' contains invalid node names. Valid node names are: {network.node_ids}.")
                if not all(state[node] in (0, 1) for node in state.keys()):
                    raise ValueError(
                        f"State '{state}' contains values other than 0 and 1.")

        else:
            raise ValueError(
                "state_space should be a list of dictionaries or a list of strings.")

        self.state_space = state_space

    def set_node_state(self, node_id, state):
        """
        Set the state of a node in all states of the state space.

        Parameters
        ----------
        node_id : str
            The identifier of the node to set the state for.
        state : int
            The state to set for the node (0 or 1).

        Raises
        ------
        ValueError
            If `node_id` is not a valid node identifier in the network, or if `state` is not 0 or 1.
        """

        if node_id not in self.network.node_ids:
            raise ValueError(
                f"Invalid node_id '{node_id}'. Valid node identifiers are: {self.network.node_ids}.")

        if state not in (0, 1):
            raise ValueError(f"Invalid state '{state}'. State must be 0 or 1.")

        for s in self.state_space:
            s[node_id] = state

    def plot_state_space(self, title="State Heatmap", cmap=None):
        """
        Plot the states of the Boolean network as a heatmap.

        Parameters
        ----------
        title : str, optional
            Title of the heatmap, default is "State Heatmap".

        cmap : ListedColormap, optional
            Colormap to use for the heatmap. Default is binary with green (1) and grey (0).

        """

        # Default binary colormap: green for 1, grey for 0
        if cmap is None:
            cmap = ListedColormap([DEFAULT_COLOURS[0], DEFAULT_COLOURS[1]])

        # Extract node names and states
        nodes = self.network.node_ids
        states = np.array([[state[node] for node in nodes] for state in self.state_space]).T

        # Create the heatmap
        plt.figure(figsize=(len(self.state_space), len(nodes)))
        plt.imshow(states, aspect="auto", cmap=cmap, interpolation="nearest")

        # Overlay "1" or "0" on the heatmap cells
        for i in range(states.shape[0]):
            for j in range(states.shape[1]):
                plt.text(j, i, str(states[i, j]), ha="center", va="center", color="black")


        # Add vertical lines between states
        for i in range(1, states.shape[1]):
            plt.axvline(i-0.5, color="white", linewidth=5)

        # remove spines around figure
        plt.gca().spines[:].set_visible(False)

        # Set axis labels
        plt.yticks(ticks=np.arange(len(nodes)), labels=[self.network.nodes[node_id].name for node_id in nodes], rotation=0)
        plt.xticks(ticks=np.arange(len(self.state_space)), labels=[f"State {i+1}" for i in range(len(self.state_space))])

        # Set title
        plt.title(title)
        plt.tight_layout()
        plt.show()

    def __repr__(self):
        return f"BooleanStateSpace(network={self.network}, state_space={self.state_space})"

class BooleanSimulationResult():
    """Class to contain simulation results"""

    def __init__(self, network, stg, initial_states):
        # TODO copy RN so that edits to original do not affect this object
        self.network = network

        self.stg = stg

        self.initial_states = initial_states

    def plot_stg(self, file=None, booldog_style=True, plot_nodes=None, use_names=True, num_characters=5):
        '''Plot the state transition graph.

        Parameters
        ----------
        fig : str
            File name of generated figure

        new_style : bool
            Whether to use booldog style, or PyBoolNet style (default) to
            plot the state transition graph. Requires pygraphviz (If not
            installed, will default to PyBoolNet version.)

        plot_nodes :  None or list of str, optional
            List of identifiers of subset of nodes to plot. If `None`, plot all nodes.
            Only valid if `booldog_style` is `True`.

        use_names : bool
            Whether to use node names instead of node identifiers in the labels.
            Only valid if `booldog_style` is `True`.


        num_characters : int
            Number of characters to use for node names in the labels.
            Only valid if `booldog_style` is `True`.

        Returns
        -------
        networkx graph or pygraphviz AGraph
            The state transition graph. If `booldog_style` is `True` and
            pygraphviz is installed, a pygraphviz AGraph object is returned.
            Otherwise, a networkx graph is returned.

        Notes
        -----
        If `booldog_style` is `True`, the nodes in the state transition graph
        are represented as tables, with each cell representing a node in the
        Boolean network. The cells are coloured green for "on" nodes and red
        for "off" nodes.
        If `use_names` is `True`, the node names are used instead of the node
        identifiers in the labels.

        '''

        if booldog_style and (importlib.util.find_spec('pygraphviz')
                              is not None):

            if plot_nodes is None:
                plot_nodes = list(self.network.node_ids)

            rev_index = {}
            labels = {}
            for node_id in self.network.node_ids:
                # build reverse index (int index -> node id)
                rev_index[self.network.index[node_id]] = node_id

                if use_names:
                    label = self.network.nodes[node_id].name
                else:
                    label = node_id
                labels[self.network.index[node_id]] = label[:num_characters]


            self.stg.graph['node']['shape'] = 'plaintext'

            for n in self.stg:
                label = '<<TABLE BORDER="1.5" CELLBORDER="1.5" CELLSPACING="2" CELLPADDING="5" COLOR="black"><TR>'
                for i, x in enumerate(n):
                    node_id = rev_index[i]
                    if node_id in plot_nodes:
                        label += f'<TD STYLE="ROUNDED" BGCOLOR="{DEFAULT_COLOURS[int(x)]}">{labels[i]}<BR/>{x}</TD>'
                label += """</TR></TABLE>>"""
                self.stg.nodes[n]["label"] = label

            agraph = nx.drawing.nx_agraph.to_agraph(self.stg)
            agraph.layout('dot')

            if file:
                agraph.draw(file)
                logger.info("Saved figure to %s", file)

            return agraph

        if booldog_style:
            logger.warning(
                "pygraphviz not available, defaulting to pyboolnet style.")

        self.stg.graph["node"]["color"] = "cyan"
        self.stg.graph["node"]["height"] = 0.3
        self.stg.graph["node"]["width"] = 0.45

        if file is not None:
            stg2image(self.stg, file, layout_engine="dot")

        return self.stg.graph



    def plot_state_space(self, title="State Heatmap", cmap=None):
        '''Plot the states of the Boolean network as a heatmap.

        Parameters
        ----------
        title : str, optional
            Title of the heatmap, default is "State Heatmap".

        cmap : ListedColormap, optional
            Colormap to use for the heatmap. Default is binary with green (1) and grey (0).

        '''

        state_space = [dict(zip(self.network.node_ids, map(int, state))) for state in self.stg.nodes()]

        BooleanStateSpace(network=self.network, state_space=state_space).plot_state_space(title=title, cmap=cmap)

    def export(self, file):
        '''Export Boolean simulation result to ??

        '''

        file_writable(file)
        #TODO

    def make_animation(self,
                 base_suid,
                 gif=None,
                 mp4=None,
                 initial_state=None,
                 colour_on=None,
                 colour_off=None,
                 cycle_repeats=3,
                 max_steps=None,
                 duration=400,
                 loop=0,
                 sleep=1):
        '''


        Parameters
        ----------

        base_suid: int
            The network SUID from Cytoscape, of the network to base the animation
            on.
        gif: None or a file path
            Path to save animation (GIF) to.
        mp4: None or a file path
            Not recommended to use, as this feature is experimental. Path to save animation (MP4) to.
        initial_state: str
            Initial state for the simulation. For valid initial states, see
            the object attributes initial_states
        colour_on: str
            Hex code for colour of "on" nodes, default="#b2df8a" (light green),
        colour_off: str
            Hex code for colour of "on" nodes, default="#6f6f6f" (light grey)
        cycle_repeats: int
            If there's a cycle in the simulation, how many times should it repeat,
            default=3,
        max_steps=None
            Maximum number of frames to include in the animation. If none, will
            be  `num_state + len(cycle)*cycle_repeats`,
            where `num_states` is the number of states in the simulation, and
            `len(cycle)` is the length of the cycle (if one exists). default=None
        duration=400
            Time on each frame (in milliseconds)
        loop: int
            Number of loops to the animation. 0 is infinite. default=0.
        sleep: int
            Seconds to sleep between Cytoscape network exports. (See notes).

        Returns
        -------


        Notes
        -----
        Do not interact with Cytoscape while the networks are being rendered, as
        this will interfere with the selection and colouring of nodes.

        Occasionally Cytoscape exports get corrupted (e.g. node borders are rendered
        in the wrong order, node fills are placed on the wrong node). This is
        independent of booldog, and it is unclear why it happens. Rerunning
        the function may help, and increasing the `sleep` parameter may also help.

        Mp4 export is in development, and not recommended to use.

        '''

        if gif is None and mp4 is None:
            raise ValueError(
                'You need to pass a file name for at least one of `gif` or `mp4`'
            )

        if not (_CYTOSCAPE_AVAILABLE and _PILLOW_AVALIBLE):
            # TODO: Exception should be split per library, leaving it for now
            raise ImportError(
                'py4cytoscape (https://py4cytoscape.readthedocs.io/) '
                'is needed to interact with Cytoscape. '
                'We suggest you install it using pip. '
                '\n\n'
                'Pillow (PIL) (https://python-pillow.org/) '
                'is needed to generate GIFs. '
                'We suggest you install it using pip. ')

        if initial_state is None:
            if len(self.initial_states) != 1:
                raise ValueError(
                    'Simulation was generated from multiple initial states. '
                    'Please specify a single initial state for the animation '
                    'using the `initial_state` argument. ')
            initial_state = self.initial_states[0]

        else:
            if not initial_state in self.initial_states:
                raise ValueError(
                    f'Passed `initial_state={initial_state}` is not a valid '
                    f'initial state for this simulation result. For valid '
                    f'initial states, see `{self.__class__}.initial_states`.')

        # state_order = dict(nx.bfs_successors(self.stg, source=initial_state))
        # all_states = [initial_state] + [n for e in list(state_order.values()) for n in e]
        state_order = {}
        state = initial_state
        while not state in state_order:
            next_states = list(self.stg.successors(state))
            if len(next_states) == 0:
                state_order[state] = None
                break
            next_state = next_states[0]
            state_order[state] = next_state
            state = next_state

        if max_steps is None:
            subgraph = nx.induced_subgraph(self.stg, state_order)
            num_states = subgraph.number_of_nodes()

            cycles = list(nx.simple_cycles(subgraph))
            if len(cycles) == 1:  # should be at most 1
                cycle_len = len(cycles[0])
            else:
                cycle_len = 0

            max_steps = num_states + cycle_repeats*cycle_len

            logger.debug(
                "Animation: num_states %s, cycle_len %s, cycle_repeats %s, max_steps %i",
                num_states, cycle_len, cycle_repeats, max_steps)

        colors = {"1": colour_on if colour_on else DEFAULT_COLOURS[1], "0": colour_off if colour_off else DEFAULT_COLOURS[0]}

        with tempfile.TemporaryDirectory() as tmpdir:

            logger.debug('Created temporary directory: %s', tmpdir)

            tmpdir_path = Path(tmpdir)

            state_images = {}

            # clones the network
            select_result = p4c.select_all_nodes(network=base_suid)
            suid = p4c.networks.create_subnetwork(
                nodes=select_result,
                subnetwork_name=f"{base_suid}-animation",
                network=base_suid)
            p4c.set_visual_style(p4c.get_current_style(base_suid),
                                 network=suid)

            for state in state_order:

                # colours the network according to the state
                for i, val in enumerate(state):
                    n = self.network.node_ids[i]
                    nodes_by_suid = p4c.select_nodes([n],
                                                     by_col="id",
                                                     network=suid)['nodes']
                    p4c.style_bypasses.set_node_color_bypass(
                        node_names=nodes_by_suid,
                        new_colors=colors[val],
                        network=suid)

                format_ = "png"
                fig_name = tmpdir_path / f"state_{state}.{format_}"
                cytoscape_utils.export_network(suid, fig_name, format=format_)

                state_images[state] = fig_name
                time.sleep(sleep)

            # TODO delete networks at suid

            current = initial_state
            images = []
            i = 0
            while i < max_steps:
                images.append(state_images[current])
                current = state_order[current]
                i += 1

            logger.debug("GIF len: %i", len(images))

            if gif:
                # use exit stack to automatically close opened images
                with contextlib.ExitStack() as stack:

                    # lazily load images
                    imgs = (stack.enter_context(Image.open(f)) for f in images)

                    # extract  first image from iterator
                    img = next(imgs)

                    # https://pillow.readthedocs.io/en/stable/handbook/image-file-formats.html#gif
                    img.save(fp=gif,
                             format='GIF',
                             append_images=imgs,
                             save_all=True,
                             duration=duration,
                             loop=loop)
                    logger.debug('Saved GIF to %s', gif)

            if mp4:

                logger.warning("Creating MP4 - this feature is experimental!")
                def animate(i):
                    im = plt.imread(images[i])
                    return [plt.imshow(im, interpolation="none")]

                fig, ax = plt.subplots()
                ax.axis("off")
                ani = animation.FuncAnimation(fig,
                                              animate,
                                              frames=len(images),
                                              interval=duration,
                                              blit=True,
                                              repeat=True,
                                              repeat_delay=duration * 4)
                ani.save(mp4)
                logger.debug('Saved MP4 to %s', mp4)
