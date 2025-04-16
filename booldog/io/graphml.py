import warnings
import logging

import xmltodict
import igraph as ig

from booldog.io.igraph import igraph2interactions

logger = logging.getLogger(__name__)

def graphml2interactions(
    file,
    edge_type_key="weight",
    node_name_key="name",
    activator_value=1,
    inhibitor_value=-1,
    yEd=False,
    inhibitor_symbol="t_shape",
    activator_symbol="standard",
    use_labels=False,
    **_
):
    '''Loads a graphml file to a dictionary of {source :{target: interaction}}.
    Since graphml is not well defined for Boolean graphs (or even standard for
    interaction networks), this import functionality has limited support.

    Parameters
    ----------
    edge_type_key: str, optional
        The edge data key that defines the interaction type as activatory or inhibitory. Only if yEd=False.
    node_name_key: str, optional
        Vertex attribute that contains the primary identifier (default="name")
    inhibitor_value: str, optional
        Value of edge_type_key for inhibition edges (default=-1). Only useful if yEd=False.
    activator_value: str, optional
        Value of edge_type_key for activation edges (default=1). Only useful if yEd=False.
    yEd: bool, optional
        If graphml file originates as a yEd export (default=False)
    inhibitor_symbol: str, optional
        Symbol of inhibition edges (default="white_diamond"). Only if yEd=True.
    activator_symbol: str, optional
        Symbol of activation edges (default="standard"). Only if yEd=True.
    '''

    with warnings.catch_warnings(record=True) as warning_list:
        g = ig.Graph.Read_GraphML(file)
    if len(warning_list) > 0:
        logger.warning(
            'Graphml parsing resulted in the following igraph warnings:\n'
            + "\n".join(["\t" + str(x.message) for x in warning_list])
        )

    if not (node_name_key in g.vertex_attributes()):
        logger.warning(
            f'"{node_name_key}" is not present as a node attribute.'
            f' Using "id" as the base name instead.'
        )
        for v in g.vs():
            v[node_name_key] = v["id"]

    if yEd:
        with open(file) as f:
            xml_dict = xmltodict.parse(f.read())

        node_labels = {}
        for v in xml_dict["graphml"]["graph"]["node"]:
            node_data = v["data"]
            if not isinstance(node_data, list):
                node_data = [node_data]
            for data_part in node_data:
                try:
                    node_labels[v["@id"]] = data_part["y:ShapeNode"]["y:NodeLabel"]["#text"]
                except KeyError:
                    pass

        if use_labels:
            # where possible, replace the node id with the node label
            for v in g.vs():
                if (v["id"] in node_labels):
                    v[node_name_key] = node_labels[v["id"]]
        else:
            # Here, note that the labels exist, and they will not be used
            # Parsing the labels is not necessary, but helpful for users not
            # aware of "id" vs "label" vs "name" conventions, especially if
            # their manual node labels go missing
            for v in g.vs():
                if (v["id"] in node_labels) and (node_labels[v["id"]] != v[node_name_key]):
                    logger.warning(
                        f"Graphml contains node labels, will not be loaded "
                        f"to the booldog.Network: Node {v['id']} has label "
                        f"{node_labels[v['id']]}."
                    )

        edge_attrs = {}
        for e in xml_dict["graphml"]["graph"]["edge"]:
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
                    f'Issue with edge "{e["@id"]}, no arrow symbol to define interaction type'
                )
            elif symbol == activator_symbol:
                edge_attrs[e["@id"]] = 1
            elif symbol == inhibitor_symbol:
                edge_attrs[e["@id"]] = -1
            else:
                logger.warning(
                    f'Issue with edge ", e["@id"] '
                    f'arrow symbol "{symbol}" not recognized as activator or '
                    f' inhibitor, perhaps you need to define the '
                    f'"inhibitor_symbol" ({inhibitor_symbol}) or '
                    f'"activator_symbol" ({activator_symbol}) in keyword '
                    f'arguments. '
                )

        for e in g.es():
            if e["id"] in edge_attrs:
                e[edge_type_key] = edge_attrs[e["id"]]
            else:
                logger.warning(f"Edge {e} has no edge type.")

    return igraph2interactions(g,
                               node_name_key=node_name_key,
                               edge_type_key=edge_type_key)
