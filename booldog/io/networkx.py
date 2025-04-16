import logging

logger = logging.getLogger(__name__)

def networkx2interactions(g,
                        edge_type_key="weight",
                        **_):
    ''' Todo
    '''
    interactions = {node: {} for node in g.nodes()}
    for source, target, data in g.edges(data=True):
        symbol = data[edge_type_key]
        interactions[source][target] = symbol

    return interactions


