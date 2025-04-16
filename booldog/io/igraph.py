import logging

logger = logging.getLogger(__name__)



def igraph2interactions(g,
                        node_name_key="name",
                        edge_type_key="weight",
                        **_):
    ''' Todo
    '''
    interactions = {node[node_name_key]: {} for node in g.vs()}
    for e in g.es():
        source = g.vs()[e.source][node_name_key]
        target = g.vs()[e.target][node_name_key]
        symbol = e[edge_type_key]

        interactions[source][target] = symbol

    return interactions
