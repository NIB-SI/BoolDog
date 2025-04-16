'''Boolean simulation results
'''

import contextlib
import importlib
import logging
import tempfile
import time
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.animation as animation
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

class BooleanSimulationResult():
    """Small class to contain simulation results"""

    def __init__(self, network, stg, initial_states):
        # TODO copy RN so that edits to original do not affect this object
        self.network = network

        self.stg = stg

        self.initial_states = initial_states

    def plot_stg(self, file=None, booldog_style=True, plot_nodes=None):
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
            Subset of nodes to plot. If `None`, plot all nodes.
            Only valid if `booldog_style` is `True`.
        '''

        if booldog_style and (importlib.util.find_spec('pygraphviz')
                              is not None):

            if plot_nodes is None:
                plot_nodes = self.network.nodes

            rev_index = {}
            for node in self.network.nodes:
                rev_index[self.network.index[node]] = node

            self.stg.graph['node']['shape'] = 'plaintext'
            colors = {"1": "#80ff8a", "0": "#ff9580"}

            for n in self.stg:
                label = '<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="1" ><TR>'
                for i, x in enumerate(n):
                    node_name = rev_index[i]
                    if node_name in plot_nodes:
                        label += f'<TD BGCOLOR="{colors[x]}">{rev_index[i][:5]} <BR/> {x} </TD>'
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

    def export(self, file):
        '''Export Boolean simulation result to ??

        '''

        file_writable(file)
        #TODO

    def make_gif(self,
                 base_suid,
                 gif=None,
                 mp4=None,
                 initial_state=None,
                 colour_on="#238B45",
                 colour_off="#CCCCCC",
                 cycle_repeats=10,
                 max_steps=None,
                 duration=400,
                 loop=5,
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
            Not usable
        initial_state: str
            Initial state for the simulation. For valid initial states, see
            the object attributes initial_states
        colour_on: str
            Hex code for colour of "on" nodes, default="#238B45" (dark green),
        colour_off: str
            Hex code for colour of "on" nodes, default="#CCCCCC" (light grey)
        cycle_repeats: int
            If there's a cycle in the simulation, how many times should it repeat,
            default=10,
        max_steps=None
            Maximum number of frames to include in the animation. If none, will
            be  `num_state + len(cycle)*cycle_repeats`,
            where `num_states` is the number of states in the simulation, and
            `len(cycle)` is the length of the cycle (if one exists). default=None
        duration=400
            Time on each frame (in milliseconds)
        loop: int
            Number of loops to the annimation. ) is infinite. default=5.
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
                "GIF num_states %s, cycle_len %s, cycle_repeats %s, max_steps %i",
                num_states, cycle_len, cycle_repeats, max_steps)

        colors = {"1": colour_on, "0": colour_off}

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
                    n = self.network.nodes[i]
                    nodes_by_suid = p4c.select_nodes([n],
                                                     by_col="name",
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
