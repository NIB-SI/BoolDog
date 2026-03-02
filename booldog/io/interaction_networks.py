''' Read standard networks formats/objects (networkx, igraph, graphml, SIF) and
convert to Boolean network.
'''

import logging
import warnings
import json
import xmltodict

try:
    import igraph as ig
    _IGRAPH_AVAILABLE = True
except ImportError as e:
    _IGRAPH_AVAILABLE = False

from booldog.classes import BoolDogNode, BoolDogModelInfo
from booldog.io.interaction_logic import interactions2rules

logger = logging.getLogger(__name__)

##############################
# Utility classes/functions
##############################


def networkx2interactions(g, edge_type_key="interaction", **_):
    ''' Convert a NetworkX Graph into a nested dictionary of interactions.

    Parameters
    ----------
    g : networkx.(Di)Graph
        The input graph object from the `Networkx` library.
        (Assumed to be directed.)

    edge_type_key : str, optional
        The edge attribute key to use for interaction values
        (e.g., weight, type). Default is "interaction".

    **_ : dict
        Additional keyword arguments (ignored).

    Returns
    -------
    dict of dict
        A nested dictionary where `interactions[source][target] = symbol`, with
        `symbol` being the edge attribute specified by `edge_type_key`.
    '''

    interactions = []
    for source, target, data in g.edges(data=True):
        symbol = data[edge_type_key]
        interactions.append((source, target, symbol))

    return interactions


def igraph2interactions(g, node_id_key="name", edge_type_key="interaction"):
    ''' Convert an igraph Graph into a list of interactions.

    Parameters
    ----------
    g : igraph.Graph
        The input graph object from the `igraph` library.
    node_id_key : str, optional
        The vertex attribute key that contains the primary identifier to use
        for node identifiers. Default is "name".
    edge_type_key : str, optional
        The edge attribute key to use for interaction type
        (e.g., weight, type). Default is "interaction".

    Returns
    -------
    dict of dict
        A nested dictionary where `interactions[source][target] = symbol`, with
        `symbol` being the edge attribute specified by `edge_type_key`.
    '''
    interactions = []
    for e in g.es():
        source = g.vs()[e.source][node_id_key]
        target = g.vs()[e.target][node_id_key]
        symbol = e[edge_type_key]

        interactions.append((source, target, symbol))

    return interactions


###############################
# In
###############################


def read_interactions(interactions_input, node_names=None, **kwargs):
    ''' Create Network from interactions (list)

    Parameters
    ----------
    interactions_input : list or str
        A list of interactions in the network. Each interaction should be a tuple of
        (source, target, sign), where source and target are node identifiers and
        sign is either activator_symbol or inhibitor_symbol.
        Or a path to a JSON file containing such a list of interactions.
    activator_symbol : int, optional
        The value representing activation in the network. Default is 1.
    inhibitor_symbol : int, optional
        The value representing inhibition in the network. Default is -1.
    **kwargs: dict
        Additional keyword arguments (ignored).

    Returns
    -------
    n: booldog.BoolDogModel
        A BoolDogModel object representing the Boolean network.

    Notes
    -----
    Uses SQUAD logic to obtain Boolean network.
    '''

    if isinstance(interactions_input, list):
        # nothing to do
        interactions = interactions_input
        source = None
    else:
        with open(interactions_input, 'r', encoding="utf-8") as f:
            interactions = json.load(f)
        source = interactions_input

    rules = interactions2rules(interactions, **kwargs)

    # all nodes (can be regulator or source)
    all_nodes = [n for e in interactions for n in e[:2]]

    nodes = [
        BoolDogNode(identifier=node_id,
                    rule=rules.get(node_id, None),
                    name=node_names.get(node_id, None) if node_names else None)
        for node_id in all_nodes
    ]

    modelinfo = BoolDogModelInfo(source=source,
                                 source_format="interactions-SQUAD")

    return {"nodes": nodes, "modelinfo": modelinfo, "primes": None}


def read_sif(file,
             delim="\t",
             source_col=0,
             target_col=1,
             interaction_col=2,
             header=True,
             **kwargs):
    '''Reads in a SIF file of interactions

    Parameters
    ----------
    file: str
        Path to the SIF file
    delim: str, optional
        Delimiter used in the SIF file (default="\t")
    header: bool, optional
        If the first line of the file is a header (default=True)
    source_col: int or str, optional
        Column index (if int) or column name (if str) of source node (default=0)
    target_col: int or str, optional
        Column index (if int) or column name (if str) of target node (default=1)
    interaction_col: int, optional
        Column index (if int) or column name (if str) of interaction type (symbol) (default=2)
    activator_symbol: str, optional
        Symbol of activation edges in `interaction_col` (default="1")
    inhibitor_symbol: str, optional
        Symbol of inhibition edges in `interaction_col` (default="-1")
    **kwargs: dict
        Additional keyword arguments (ignored).

    Returns
    -------


    Notes
    -----
    Uses SQUAD logic to obtain Boolean network.
    '''

    interactions = []
    with open(file, "r", encoding="utf-8") as handle:
        if header:
            columns = handle.readline().rstrip().split(delim)
            if isinstance(source_col, str):
                source_col = columns.index(source_col)
            if isinstance(target_col, str):
                target_col = columns.index(target_col)
            if isinstance(interaction_col, str):
                interaction_col = columns.index(interaction_col)

        for line in handle:
            cols = line.rstrip().split(delim)
            source, target, interaction = cols[source_col], cols[
                target_col], cols[interaction_col]
            interactions.append((source, target, interaction))

    rules = interactions2rules(interactions, **kwargs)

    # all nodes (can be regulator or source)
    all_nodes = [n for e in interactions for n in e[:2]]

    nodes = [
        BoolDogNode(
            identifier=node_id,
            rule=rules.get(node_id, None),
            # name=node_names.get(node_id, None) if node_names else None
        ) for node_id in all_nodes
    ]

    modelinfo = BoolDogModelInfo(source=file,
                                 source_format="interactions-SQUAD")

    return {"nodes": nodes, "modelinfo": modelinfo, "primes": None}


def read_igraph(g,
                node_id_key="name",
                edge_type_key="interaction",
                node_name_key="name",
                **kwargs):
    ''' Create BooleanNetwork from a igraph Graph object.

    Parameters
    ----------
    g: igraph.Graph
        The input graph object from the `igraph` library.
    node_id_key : str, optional
        The vertex attribute key that contains the primary identifier to use for node names. Default is "name".
    node_name_key : str, optional
        The vertex attribute key that contains the node name (e.g. display label). Default is "name".
    edge_type_key : str, optional
        The edge attribute key to use for interaction values (e.g., weight, type). Default is "interaction".
    activator_symbol: str, optional
        Symbol or value of activation edges in `edge_type_key` of g (default="1")
    inhibitor_symbol: str, optional
        Symbol or value of inhibition edges in `edge_type_key` of g (default="-1")

    **kwargs: dict
        Additional keyword arguments (ignored).

    Returns
    -------



    Notes
    -----
    Uses SQUAD logic to obtain Boolean network.
    '''

    if not _IGRAPH_AVAILABLE:
        raise ImportError("igraph is not available.")

    interactions = igraph2interactions(g,
                                       node_id_key=node_id_key,
                                       edge_type_key=edge_type_key)

    rules = interactions2rules(interactions, **kwargs)

    nodes = [
        BoolDogNode(identifier=node[node_id_key],
                    rule=rules.get(node[node_id_key], None),
                    name=node[node_name_key]
                    if node_name_key in g.vertex_attributes() else None)
        for node in g.vs()
    ]

    modelinfo = BoolDogModelInfo(source="object", source_format="igraph-SQUAD")

    return {"nodes": nodes, "modelinfo": modelinfo, "primes": None}


def read_networkx(g,
                  edge_type_key="interaction",
                  node_name_key="name",
                  **kwargs):
    ''' Create BooleanNetwork from a networkx (Di)Graph object.

    Parameters
    ----------
    g: networkx.(Di)Graph
        The input graph object from the `networkx` library.
    edge_type_key : str, optional
        The edge attribute key to use for interaction values (e.g., weight, type). Default is "interaction".
    activator_symbol: str, optional
        Symbol or value of activation edges in `edge_type_key` of g (default="1")
    inhibitor_symbol: str, optional
        Symbol or value of inhibition edges in `edge_type_key` of g (default="-1")
    node_name_key : str, optional
        The vertex attribute key that contains the node name (display label). Default is "name".

    **kwargs: dict
        Additional keyword arguments (ignored).

    Returns
    -------
    n: booldog.BoolDogModel
        A BoolDogModel object representing the Boolean network.

    Notes
    -----
    Uses SQUAD logic to obtain Boolean network.
    '''
    interactions = networkx2interactions(g, **kwargs)

    rules = interactions2rules(interactions, **kwargs)

    nodes = [
        BoolDogNode(identifier=node_id,
                    rule=rules.get(node_id, None),
                    name=data.get("name", None))
        for node_id, data in g.nodes(data=True)
    ]

    modelinfo = BoolDogModelInfo(source="object",
                                 source_format="networkx-SQUAD")

    return {"nodes": nodes, "modelinfo": modelinfo, "primes": None}


def read_graphml(file,
                 edge_type_key="interaction",
                 node_id_key="name",
                 node_name_key="name",
                 yEd_labels=False,
                 yEd_arrow_head=False,
                 use_labels=True,
                 **kwargs):
    ''' Extract relavent parts for a Boolean network from a graphml file.

    Since graphml is not well defined for Boolean networks (or even standard for
    interaction networks), this read functionality has limited support.

    Parameters
    ----------
    file: str
        Path to the graphml file
    edge_type_key: str, optional
        The edge attribute key to use for interaction type (e.g., weight, type). Default is "interaction".
        Only if yEd=False (yEd uses the arrow head symbol).
    node_id_key : str, optional
        The vertex attribute key that contains the primary identifier to use for node names. Default is "name".
    node_name_key : str, optional
        The vertex attribute key that contains the node name (e.g. display label). Default is "name".
    yEd_labels: bool, optional
        If graphml file originates as a yEd export, the node attribute
        "y:NodeLabel" is used to determine node names (default=False)
    yEd_arrow_head: bool, optional
        If graphml file originates as a yEd export, the edge attribute
        "y:Arrows" is used to determine interaction type (default=False)
    activator_symbol: int or str, optional
        Symbol or value of activation edges. If yEd=True, default="standard", else default is `1`.
    inhibitor_symbol: int or str, optional
        Symbol or value of inhibition edges. If yEd=True, default="t_shape", else default is `-1`.

    Returns
    -------
    n: booldog.BoolDogModel
        A BoolDogModel object representing the Boolean network.

    Notes
    -----
    Uses igraph to parse the graphml file, and extracts node and edge
    attributes to determine interactions.

    If `yEd_label=True`, it also parses the yEd-specific attributes to
    extract node labels (y:NodeLabel) for node names.

    If `yEd_arrow_head=True`, it parses the yEd-specific edge arrow types
    (arrow head symbols of y:Arrows) to determine interaction types
    (activation/inhibition).
    '''

    if not _IGRAPH_AVAILABLE:
        raise ImportError("igraph is not available.")

    # decide on activation/inhibition symbols
    if yEd_arrow_head:
        activator_symbol = kwargs.get("activator_symbol", "standard")
        inhibitor_symbol = kwargs.get("inhibitor_symbol", "t_shape")
    else:
        activator_symbol = kwargs.get("activator_symbol", 1)
        inhibitor_symbol = kwargs.get("inhibitor_symbol", -1)

    with warnings.catch_warnings(record=True) as warning_list:
        g = ig.Graph.Read_GraphML(file)

        if len(warning_list) > 0:
            logger.warning(
                'Graphml parsing resulted in the following igraph warnings:\n%s',
                "\n".join(["\t" + str(x.message) for x in warning_list]))

    if yEd_labels or yEd_arrow_head:
        with open(file, encoding="utf-8") as f:
            xml_dict = xmltodict.parse(f.read())

        if yEd_labels:
            for v in xml_dict["graphml"]["graph"]["node"]:
                graphml_node_id = v["@id"]
                label = None
                node_data = v["data"]
                if not isinstance(node_data, list):
                    node_data = [node_data]
                ynode = None
                for data_part in node_data:
                    if "y:ShapeNode" in data_part:
                        ynode = data_part['y:ShapeNode']
                        break
                if ynode and "y:NodeLabel" in ynode.keys():
                    if "#text" in ynode["y:NodeLabel"]:
                        label = ynode["y:NodeLabel"]["#text"]
                else:
                    logger.debug("No y:NodeLabel for node %s", graphml_node_id)

                if label is not None:
                    g.vs.find(id=graphml_node_id)[node_name_key] = label
                else:
                    logger.warning("Node %s has no y:NodeLabel.",
                                   graphml_node_id)

        if yEd_arrow_head:
            edge_attrs = {}
            for e in xml_dict["graphml"]["graph"]["edge"]:
                graphml_edge_id = e["@id"]
                symbol = None
                edge_data = e["data"]
                if not isinstance(edge_data, list):
                    edge_data = [edge_data]
                for data_part in edge_data:
                    try:
                        symbol = data_part["y:PolyLineEdge"]['y:Arrows']['@target']
                    except KeyError:
                        pass
                if symbol is None:
                    logger.warning(
                        'Issue with Graphml edge id "%s", no arrow symbol to define interaction type',
                        graphml_edge_id)

                logger.debug("Edge %s has yEd arrow type: %s", e["@id"],
                             symbol)

                # This does not work with igraph=1.0.0 since it does
                # not correctly import edge ids:
                # https://github.com/igraph/igraph/issues/2892
                # Updated once issue fixed, to eliminate the for loop below
                # if symbol is not None:
                #     g.es.find(id=graphml_edge_id)[edge_type_key] = symbol

                # instead, to prevent Error
                edge_attrs[graphml_edge_id] = symbol

            for e in g.es():
                if e["id"] in edge_attrs:
                    e[edge_type_key] = edge_attrs[e["id"]]
                else:
                    e[edge_type_key] = None
                    logger.warning("Edge %s has no edge type.", e["id"])

    # check if node_id_key is present
    if not (node_id_key in g.vertex_attributes()):
        logger.warning(
            '"%s" is not present as a node attribute.'
            ' Using "id" as the primary identifier instead.', node_id_key)
        for v in g.vs():
            v[node_id_key] = v["id"]

    # check if node_name_key is present
    # if its missing, fall back to node_id_key, log to user
    if node_name_key != node_id_key and not (node_name_key
                                             in g.vertex_attributes()):
        logger.warning(
            '"%s" is not present as a node attribute. Node names will be the same as node identifiers.',
            node_name_key)
        for v in g.vs():
            v[node_name_key] = v[node_id_key]

    interactions = igraph2interactions(g,
                                       node_id_key=node_id_key,
                                       edge_type_key=edge_type_key)

    rules = interactions2rules(interactions,
                               activator_symbol=activator_symbol,
                               inhibitor_symbol=inhibitor_symbol)

    nodes = [
        BoolDogNode(identifier=node[node_id_key],
                    rule=rules.get(node[node_id_key], None),
                    name=node[node_name_key]
                    if node_name_key in g.vertex_attributes() else None)
        for node in g.vs()
    ]

    modelinfo = BoolDogModelInfo(source=file, source_format="graphml-SQUAD")

    return {"nodes": nodes, "modelinfo": modelinfo, "primes": None}
