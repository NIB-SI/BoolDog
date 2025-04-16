'''Booldog networks interact with Cytoscape
'''


try:
    import py4cytoscape as p4c
    from booldog.utils import cytoscape_utils
    _CYTOSCAPE_AVAILABLE = True
except ImportError as e:
    _CYTOSCAPE_AVAILABLE = False

from booldog.resources import cytoscape_style_xml



def booldog2cytoscape(network, as_logic_circuit=False,
    collection="Booldog Network", layout=None, style=None
    ):

    g = network.to_networkx(as_logic_circuit=as_logic_circuit)

    suid = p4c.networks.create_network_from_networkx(g, collection=collection)

    if layout:
        p4c.layouts.layout_network(layout_name=layout, network=suid)

    if not style:
        style = "logical_circuit"
        if not style in p4c.styles.get_visual_style_names():
            p4c.import_visual_styles(cytoscape_style_xml)

    p4c.styles.set_visual_style(style, network=suid)

    return suid