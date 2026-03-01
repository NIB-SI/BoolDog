'''Booldog networks interact with Cytoscape
'''

from booldog.resources import cytoscape_style_xml
from booldog.utils.decorators import silence_stdout

try:
    import py4cytoscape as p4c
    from booldog.utils import cytoscape_utils
    _CYTOSCAPE_AVAILABLE = True
except ImportError as e:
    _CYTOSCAPE_AVAILABLE = False

# ---- Utility functions ------------------------------------------------------

def silence_p4c_loggers():
    import logging
    loggers = [
        logging.getLogger("py4..."),
        logging.getLogger("py4...S"),
    ]
    for logger in loggers:
        logger.setLevel(logging.CRITICAL)

if _CYTOSCAPE_AVAILABLE:
    silence_p4c_loggers()

@silence_stdout
def test_cytoscape_connection():
    """Test if a connection to Cytoscape can be established."""
    if not _CYTOSCAPE_AVAILABLE:
        raise ImportError(
            'py4cytoscape (https://py4cytoscape.readthedocs.io/) '
            'is needed to interact with Cytoscape. '
            'We suggest you install it using pip. ')
    try:
        p4c.cytoscape_ping()
        return True
    except Exception as e:
        raise ConnectionError(
            "Unable to connect to Cytoscape. Make sure Cytoscape is running.") from e


# ---- Main function ----------------------------------------------------------

def booldog2cytoscape(model, as_logic_circuit=False, title=None,
    collection="Booldog Network", layout=None, style=None
    ):
    """Display a BoolDog Boolean model in Cytoscape.

    Parameters
    ----------
    model : booldog:BoolDogModel
        A BoolDog object representing a Boolean network.

    as_logic_circuit: bool
        If True, the graph is exported as a logic circuit (Boolean rules
        are represented as "logical" nodes (and, or, not) and edges.
        Otherwise, it is exported as a interaction graph. Default is False.

    title: str or None
        The title of the Cytoscape network. Default is '{model_id} interaction network' if
        `as_logic_circuit` is False and '{model_id} logic circuit' else.

    collection: str
        The name of the Cytoscape collection to add the network to. Default is "Booldog Network".

    layout: str or None
        The name of the Cytoscape layout to apply to the network. Default is None (no layout).

    style: str or None
        The name of the Cytoscape visual style to apply to the network.
        If None, a default style will be applied. Default is None.

    Returns
    -------
    suid : int
        The SUID of the created Cytoscape network.
    """

    test_cytoscape_connection()

    if title is None:
        title = f"{model.modelinfo.identifier} logic circuit" if as_logic_circuit else f"{model.modelinfo.identifier} interaction network"

    g = model.to_networkx(as_logic_circuit=as_logic_circuit)

    # pyboolnet edges have an attribute "sign" which is a set of 1 (activation)
    # or -1 (inhibition). These are not supported by Cytoscape, so
    # overwrite to a str (and use the "type" edge attribute)
    if not as_logic_circuit:
        for _, _, data in g.edges(data=True):
            data["sign"] = str(data['sign'])

    # add a label attribute to each node for Cytoscape style
    for node_id in model.node_ids:
        g.nodes[node_id]["label"] = model.nodes[node_id].name

    suid = p4c.networks.create_network_from_networkx(g, title=title, collection=collection)

    if layout:
        p4c.layouts.layout_network(layout_name=layout, network=suid)

    if not style:
        style = "logical_circuit"
        if not style in p4c.styles.get_visual_style_names():
            p4c.import_visual_styles(cytoscape_style_xml)

        p4c.styles.set_visual_style(style, network=suid)

    return suid
