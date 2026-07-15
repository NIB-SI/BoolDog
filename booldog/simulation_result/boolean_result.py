'''Boolean simulation results

Classes returned by :py:meth:`~booldog.boolean.boolean.BooleanNetworkMixin.boolean_simulation`
(:class:`BooleanSimulationResult`, wrapping the state transition graph) and used to
represent/plot subsets of the Boolean state space (:class:`BooleanStateSpace`).
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
'''dict of int -> str : Default hex colours used to represent node states of
1 ("on", light green) and 0 ("off", medium grey) in heatmaps, state
transition graph tables, and animations.'''

class BooleanStateSpace():
    """Class representing a subspace (a collection of states) of a Boolean
    network, e.g. the state space explored by a
    :class:`BooleanSimulationResult`, or an arbitrary user-supplied set of
    states."""

    def __init__(self, network, state_space):
        """
        Initialize the BooleanState object.

        Parameters
        ----------
        network : BoolDogModel
            The Boolean network object.
        state_space : list of dict or list of str
            Either: A list of dictionaries representing the state space of the network.
            Each dictionary contains node names as keys and their states (0 to n) as values.
            Or: A list of strings representing the state space of the network,
            where each string is a binary representation of the state
            (e.g. "1010" for a network with 4 nodes).
        """
        self.network = network
        '''BoolDogModel : the Boolean network this state space belongs to.'''

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
        '''list of dict : the states of the state space, each a mapping of
        node identifier to state (0 or 1).'''

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

    def plot_state_space(self, title="State Heatmap", plot_nodes=None, cmap=None):
        """
        Plot the states of the Boolean network as a heatmap.

        Parameters
        ----------
        title : str, optional
            Title of the heatmap, default is "State Heatmap".
        plot_nodes: None or list of str, optional
            Subset of nodes to plot. If `None`, plot all nodes.
        cmap : ListedColormap, optional
            Colormap to use for the heatmap. Default is binary with green (1) and grey (0).

        Returns
        -------
        None
            Displays the heatmap directly via ``plt.show()``; the figure is
            not returned or saved to file.

        """

        # Default binary colormap: green for 1, grey for 0
        if cmap is None:
            cmap = ListedColormap([DEFAULT_COLOURS[0], DEFAULT_COLOURS[1]])

        # Extract node names and states
        if plot_nodes:
            nodes = [node_id for node_id in plot_nodes if node_id in self.network.node_ids]
        else:
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
        '''str : ``BooleanStateSpace(network=..., state_space=...)``.'''
        return f"BooleanStateSpace(network={self.network}, state_space={self.state_space})"

class BooleanSimulationResult():
    """Class to contain the result of a Boolean (synchronous) simulation:
    the explored state transition graph together with the initial state(s)
    it was generated from. Returned by
    :py:meth:`~booldog.boolean.boolean.BooleanNetworkMixin.boolean_simulation`."""

    def __init__(self, network, stg, initial_states):
        # TODO copy RN so that edits to original do not affect this object
        '''
        Parameters
        ----------
        network : BoolDogModel
            The Boolean network the simulation was run on.
        stg : networkx.DiGraph
            The state transition graph, as generated by
            ``pyboolnet.state_transition_graphs.primes2stg``. Node labels are
            strings of ``'0'``/``'1'`` (one character per node, in the order
            of ``network.node_ids``).
        initial_states : list of str or None
            The initial state(s) the state transition graph was generated
            from, in the same string encoding as the ``stg`` node labels.
        '''
        self.network = network
        '''BoolDogModel : the Boolean network the simulation was run on.'''

        self.stg = stg
        '''networkx.DiGraph : the state transition graph, with node labels
        as ``'0'``/``'1'`` strings (see :class:`BooleanSimulationResult`).'''

        self.initial_states = initial_states
        '''list of str or None : the initial state(s) used to generate the
        state transition graph.'''

    def plot_stg(self, file=None, booldog_style=True, plot_nodes=None, use_names=True, num_characters=5):
        '''Plot the state transition graph.

        Parameters
        ----------
        file : str or None, optional
            File name to save the generated figure to. If `None` (default),
            the figure is not saved to file.

        booldog_style : bool, optional
            Whether to use booldog style (default, `True`) or PyBoolNet style
            to plot the state transition graph. The booldog style requires
            pygraphviz; if it is not installed, falls back to the PyBoolNet
            style with a warning.

        plot_nodes :  None or list of str, optional
            List of identifiers of subset of nodes to plot. If `None`, plot all nodes.
            Only valid if `booldog_style` is `True`.

        use_names : bool, optional
            Whether to use node names instead of node identifiers in the labels.
            Only valid if `booldog_style` is `True`. Default `True`.

        num_characters : int, optional
            Number of characters to truncate node names/identifiers to in the
            labels. Only valid if `booldog_style` is `True`. Default 5.

        Returns
        -------
        pygraphviz.AGraph or dict
            If `booldog_style` is `True` and pygraphviz is available, the
            pygraphviz ``AGraph`` built from the state transition graph
            (with table-styled node labels) is returned. Otherwise, the
            ``graph`` attribute dictionary of the underlying networkx state
            transition graph (``self.stg.graph``) is returned.

        Notes
        -----
        If `booldog_style` is `True`, the nodes in the state transition graph
        are represented as tables, with each cell representing a node in the
        Boolean network. The cells are coloured using :data:`DEFAULT_COLOURS`
        (green for "on"/1, grey for "off"/0).
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
        '''Plot the states visited in the state transition graph as a heatmap.

        Builds a :class:`BooleanStateSpace` from every node (state) currently
        in ``self.stg`` and delegates to its
        :py:meth:`BooleanStateSpace.plot_state_space`.

        Parameters
        ----------
        title : str, optional
            Title of the heatmap, default is "State Heatmap".

        cmap : ListedColormap, optional
            Colormap to use for the heatmap. Default is binary with green (1) and grey (0).

        Returns
        -------
        None
            Displays the heatmap directly via ``plt.show()``.
        '''

        state_space = [dict(zip(self.network.node_ids, map(int, state))) for state in self.stg.nodes()]

        BooleanStateSpace(network=self.network, state_space=state_space).plot_state_space(title=title, cmap=cmap)

    def export(self, file):
        '''Export the Boolean simulation result to a file.

        Notes
        -----
        Not yet implemented: this currently does nothing.

        Parameters
        ----------
        file : str or Path
            Path to the output file.
        '''
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
        '''Render an animated GIF and/or MP4 of the trajectory from
        `initial_state`, coloured over a live Cytoscape network view.

        Requires ``py4cytoscape`` (a running Cytoscape session with the
        network already loaded) and ``Pillow``.

        Parameters
        ----------
        base_suid : int
            The network SUID from Cytoscape, of the network to base the animation
            on.
        gif : str or Path or None, optional
            Path to save animation (GIF) to. At least one of `gif`/`mp4` must
            be given.
        mp4 : str or Path or None, optional
            Not recommended to use, as this feature is experimental. Path to
            save animation (MP4) to. At least one of `gif`/`mp4` must be given.
        initial_state : str or None, optional
            Initial state to start the animated trajectory from, as a
            ``'0'``/``'1'`` string matching `self.stg` node labels. For valid
            initial states, see the object attribute `initial_states`. If
            `None` (default), uses the single state in `self.initial_states`
            (raises `ValueError` if there is more than one).
        colour_on : str or None, optional
            Hex code for colour of "on" (1) nodes, default `None` uses
            ``"#b2df8a"`` (light green).
        colour_off : str or None, optional
            Hex code for colour of "off" (0) nodes, default `None` uses
            ``"#6f6f6f"`` (medium grey).
        cycle_repeats : int, optional
            If there's a cycle in the trajectory, how many times it should
            repeat when computing the default `max_steps`. Only used if
            `max_steps` is `None`. Default 3.
        max_steps : int or None, optional
            Maximum number of frames to include in the animation. If `None`
            (default), computed as ``num_states + cycle_len * cycle_repeats``,
            where `num_states` is the number of distinct states reached from
            `initial_state`, and `cycle_len` is the length of the cycle
            reached (0 if the trajectory does not cycle back on itself).
        duration : int, optional
            Time on each frame, in milliseconds. Default 400. Also used to
            derive the `repeat_delay` for the MP4 animation (``duration * 4``).
        loop : int, optional
            Number of times the GIF should loop. 0 is infinite. Default 0.
            Only applies to `gif` output.
        sleep : int, optional
            Seconds to sleep between Cytoscape network image exports (one per
            state visited). Default 1. (See Notes.)

        Returns
        -------
        None
            The animation is written directly to `gif` and/or `mp4`; nothing
            is returned.

        Raises
        ------
        ValueError
            If neither `gif` nor `mp4` is given, or if `initial_state` is
            given but is not a member of `self.initial_states`, or if
            `initial_state` is `None` and `self.initial_states` does not
            contain exactly one state.
        ImportError
            If ``py4cytoscape`` and/or ``Pillow`` are not installed.

        Notes
        -----
        The animated trajectory is built by repeatedly following the first
        successor of each state in `self.stg` starting from `initial_state`;
        it is intended for deterministic (e.g. synchronous) trajectories.

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
